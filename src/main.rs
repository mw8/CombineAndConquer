use bevy::prelude::*;
use bevy_prototype_lyon::prelude::*;
use rand::prelude::*;
use voronoice::*;

fn main() {
    App::new()
        .insert_resource(Msaa { samples: 4 })
        .add_plugins(DefaultPlugins)
        .add_plugin(ShapePlugin)
        .add_startup_system(setup_system)
        .run();
}

// push element x onto vector v if x is already not in v
fn push_if_unique(v: &mut Vec<usize>, x: usize) {
    if !v.iter().any(|y| *y == x) {
        v.push(x);
    }
}

// map a number between -1 and 1 to an interval of length l centered at zero
fn map_interval(x: f64, l: f64) -> f64 {
    if x >= 0.0 {
        x.powf(0.75) * l * 0.5
    } else {
        -(-x).powf(0.75) * l * 0.5
    }
}

#[derive(Clone, Copy)]
enum CellType {
    Empty,         // ignore this cell
    Precinct(i32), // precinct with given vote margin
    Land,          // non-precinct land cell
    Sea,           // non-precinct sea cell
}

struct GameMap {
    rng: ThreadRng,
    width: f64,
    height: f64,
    sites: Vec<Vec2>,           // point in the middle of each cell
    types: Vec<CellType>,       // type of each cell
    vertices: Vec<Vec<Vec2>>,   // vertices of each cell
    neighbors: Vec<Vec<usize>>, // neighbors of each cell
    merged: Vec<bool>,          // whether a given cell has been merged
    districts: Vec<Vec<usize>>, // collections of cells to form districts
}

impl GameMap {
    // create a new random game map
    fn new(num_cells: usize, width: f64, height: f64) -> GameMap {
        // random starting sites for cells
        let mut rng = thread_rng();
        // let mut rng = ChaCha8Rng::seed_from_u64(0);
        let mut rand_sites = Vec::<Point>::with_capacity(num_cells);
        for _ in 0..num_cells {
            let x = map_interval(rng.gen_range(-1.0..1.0), width);
            let y = map_interval(rng.gen_range(-1.0..1.0), height);
            rand_sites.push(Point { x: x, y: y })
        }
        // build voronoi cells from the set of starting sites
        let voronoi = VoronoiBuilder::default()
            .set_sites(rand_sites)
            .set_bounding_box(BoundingBox::new_centered(width, height))
            .set_lloyd_relaxation_iterations(3)
            .build()
            .unwrap();
        // extract voronoi cell properties
        let mut sites: Vec<Vec2> = Vec::with_capacity(num_cells);
        let mut types: Vec<CellType> = Vec::with_capacity(num_cells);
        let mut vertices: Vec<Vec<Vec2>> = Vec::with_capacity(num_cells);
        for cell in voronoi.iter_cells() {
            let mut v: Vec<Vec2> = Vec::with_capacity(cell.iter_vertices().count());
            for vertex in cell.iter_vertices() {
                v.push(Vec2::new(vertex.x as f32, vertex.y as f32));
            }
            vertices.push(v);
            let cell_site = cell.site_position();
            sites.push(Vec2::new(cell_site.x as f32, cell_site.y as f32));
            if cell.is_on_hull() {
                types.push(CellType::Land);
            } else {
                types.push(CellType::Precinct(0));
            }
        }
        assert!(sites.len() == num_cells);
        assert!(types.len() == num_cells);
        assert!(vertices.len() == num_cells);
        // extract neighbor information using associated Delaunay triangulation
        let mut neighbors: Vec<Vec<usize>> = Vec::with_capacity(num_cells);
        for _ in 0..num_cells {
            neighbors.push(Vec::new());
        }
        let num_triangles = voronoi.triangulation().len();
        for i in 0..num_triangles {
            let i0 = voronoi.triangulation().triangles[3 * i + 0];
            let i1 = voronoi.triangulation().triangles[3 * i + 1];
            let i2 = voronoi.triangulation().triangles[3 * i + 2];
            push_if_unique(&mut neighbors[i0], i1);
            push_if_unique(&mut neighbors[i0], i2);
            push_if_unique(&mut neighbors[i1], i0);
            push_if_unique(&mut neighbors[i1], i2);
            push_if_unique(&mut neighbors[i2], i0);
            push_if_unique(&mut neighbors[i2], i1);
        }
        assert!(neighbors.len() == num_cells);
        GameMap {
            rng: rng,
            width: width,
            height: height,
            sites: sites,
            types: types,
            vertices: vertices,
            neighbors: neighbors,
            merged: vec![false; num_cells],
            districts: Vec::new(),
        }
    }

    // compute number of cells
    fn num_cells(&self) -> usize {
        self.types.len()
    }
    // compute number of cells that are precincts
    fn num_precincts(&self) -> usize {
        self.types
            .iter()
            .map(|&x| matches!(x, CellType::Precinct(_)) as usize)
            .sum::<usize>()
    }

    // output vector of precinct vertices
    fn precinct_vertices(&self) -> Vec<Vec<Vec2>> {
        let num_cells = self.num_cells();
        let mut vertices = Vec::new();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Precinct(_) => vertices.push(self.vertices[i].clone()),
                _ => {}
            }
        }
        vertices
    }

    // output vote margin numbers for each precinct
    fn precinct_margins(&self) -> Vec<i32> {
        let num_cells = self.num_cells();
        let mut margins = Vec::new();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Precinct(x) => margins.push(x),
                _ => {}
            }
        }
        margins
    }

    // output site location (barycenter) for each precinct
    fn _precinct_sites(&self) -> Vec<Vec2> {
        let num_cells = self.num_cells();
        let mut sites = Vec::new();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Precinct(_) => sites.push(self.sites[i]),
                _ => {}
            }
        }
        sites
    }

    // output cell index for each precinct
    fn _precinct_indices(&self) -> Vec<usize> {
        let num_cells = self.num_cells();
        let mut indices = Vec::new();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Precinct(_) => indices.push(i),
                _ => {}
            }
        }
        indices
    }

    // output vector of non-precinct land vertices
    fn land_vertices(&self) -> Vec<Vec<Vec2>> {
        let num_cells = self.num_cells();
        let mut vertices = Vec::new();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Land => vertices.push(self.vertices[i].clone()),
                _ => {}
            }
        }
        vertices
    }

    // output vector of non-precinct sea vertices
    fn sea_vertices(&self) -> Vec<Vec<Vec2>> {
        let num_cells = self.num_cells();
        let mut vertices = Vec::new();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Sea => vertices.push(self.vertices[i].clone()),
                _ => {}
            }
        }
        vertices
    }

    // remove a random precinct near the edge
    fn remove_precinct_edge(&mut self) {
        let num_cells = self.num_cells();
        let mut i0 = self.rng.gen_range(0..num_cells);
        let mut j0 = 0;
        while j0 < num_cells {
            if matches!(self.types[i0], CellType::Land) {
                let neighbors = &self.neighbors[i0];
                let num_neighbors = neighbors.len();
                let mut i1 = self.rng.gen_range(0..num_neighbors);
                let mut j1 = 0;
                while j1 < num_neighbors {
                    if matches!(self.types[neighbors[i1]], CellType::Precinct(_)) {
                        self.types[neighbors[i1]] = CellType::Land;
                        return;
                    }
                    i1 = (i1 + 1) % num_neighbors;
                    j1 += 1;
                }
            }
            i0 = (i0 + 1) % num_cells;
            j0 += 1;
        }
    }

    // change all connected precinct components smaller than a threshold to not be precincts
    fn remove_precinct_comp(&mut self, threshold: usize) {
        let num_cells = self.num_cells();
        let mut cell_checked = vec![false; num_cells];
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Empty | CellType::Land | CellType::Sea => {
                    cell_checked[i] = true;
                }
                _ => {}
            }
        }
        while cell_checked.iter().any(|x| !*x) {
            // find index of first unchecked cell
            let first_idx = cell_checked.iter().position(|&x| !x).unwrap();
            // compute size of component
            let mut in_component = vec![false; num_cells];
            in_component[first_idx] = true;
            cell_checked[first_idx] = true;
            let mut buf = vec![first_idx];
            while let Some(idx) = buf.pop() {
                for &k in &self.neighbors[idx] {
                    if !in_component[k] && matches!(self.types[k], CellType::Precinct(_)) {
                        in_component[k] = true;
                        cell_checked[k] = true;
                        buf.push(k);
                    }
                }
            }
            let component_size = in_component.iter().map(|&x| x as i32).sum::<i32>() as usize;
            if component_size < threshold {
                for i in 0..num_cells {
                    if in_component[i] {
                        self.types[i] = CellType::Land;
                    }
                }
            }
        }
    }

    fn merge_cells(&mut self, i0: usize, i1: usize) -> Result<(), &'static str> {
        let mut v0 = self.vertices[i0].clone();
        let mut v1 = self.vertices[i1].clone();
        // mark second cell as empty
        self.types[i1] = CellType::Empty;
        // find first vertex v in common between two cells
        let mut v = Vec2 {
            x: std::f32::NAN,
            y: std::f32::NAN,
        };
        let mut found = false;
        for &a0 in &v0 {
            for &a1 in &v1 {
                if a0 == a1 {
                    v = a0;
                    found = true;
                    break;
                }
            }
            if found {
                break;
            }
        }
        if !found {
            return Err("Cannot merge cells because they have no vertices in common.");
        }
        let mut v_idx0 = v0.iter().position(|&x| x == v).unwrap();
        let mut v_idx1 = v1.iter().position(|&x| x == v).unwrap();
        // let w be successor of v in v0 and the predecessor of v in v1
        let mut w_idx0 = (v_idx0 + v0.len() + 1) % v0.len();
        let w_idx1 = (v_idx1 + v1.len() - 1) % v1.len();
        // if w does not match between v0 and v1, then set w to v and set v to its predecessor
        if v0[w_idx0] != v1[w_idx1] {
            w_idx0 = v_idx0;
            // w_idx1 = v_idx1;
            v_idx0 = (v_idx0 + v0.len() - 1) % v0.len();
            v_idx1 = (v_idx1 + v1.len() + 1) % v1.len();
            v = v0[v_idx0];
            assert!(v == v1[v_idx1]);
        }
        // merge polygons
        let mut v01 = v0.split_off(w_idx0);
        let mut v11 = v1.split_off(v_idx1);
        if !v0.is_empty() {
            v0.pop();
        } else {
            v01.pop();
        }
        if !v1.is_empty() {
            v1.pop();
        } else {
            v11.pop();
        }
        v0.append(&mut v11);
        v0.append(&mut v1);
        v0.append(&mut v01);
        self.vertices[i0] = v0;
        self.vertices[i1].clear();
        self.merged[i0] = true;
        self.merged[i1] = true;
        self.sites[i0] = (self.sites[i0] + self.sites[i1]) * 0.5;
        Ok(())
    }

    // combine two random neighboring precincts
    fn combine_precincts(&mut self) -> Result<(), &'static str> {
        // find random precinct
        let num_cells = self.num_cells();
        let mut i0 = self.rng.gen_range(0..num_cells);
        let mut j = 0;
        while j < num_cells && (self.merged[i0] || !matches!(self.types[i0], CellType::Precinct(_)))
        {
            i0 = (i0 + 1) % num_cells;
            j += 1;
        }
        if j == num_cells {
            return Err("Unable to find any precincts to combine.");
        }
        // find a neighboring precinct
        let neighbors = &self.neighbors[i0];
        let num_neighbors = neighbors.len();
        let mut i1 = self.rng.gen_range(0..num_neighbors);
        j = 0;
        while j < num_neighbors
            && (self.merged[neighbors[i1]]
                || !matches!(self.types[neighbors[i1]], CellType::Precinct(_)))
        {
            i1 = (i1 + 1) % num_neighbors;
            j += 1;
        }
        if j == num_neighbors {
            return Err("Unable to find any neighbors of the selected precinct.");
        }
        i1 = neighbors[i1];
        // merge the precincts
        self.merge_cells(i0, i1)?;
        // success
        Ok(())
    }

    // randomly assign vote margins to each precinct
    fn randomize_precincts(&mut self) {
        let perlin = perlin_noise::PerlinNoise::new();
        let mut perlin_vals = Vec::new();
        let mut perlin_idxs = Vec::new();
        let mut j = 0;
        let num_cells = self.types.len();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Precinct(_) => {
                    perlin_vals
                        .push(perlin.get2d([self.sites[i].x as f64, self.sites[i].y as f64]));
                    perlin_idxs.push(j);
                    j += 1;
                }
                _ => {}
            }
        }
        perlin_idxs.sort_by(|&a, &b| perlin_vals[a].partial_cmp(&perlin_vals[b]).unwrap());
        let idx_cutoff = perlin_vals.len() / 2;
        j = 0;
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Precinct(_) => {
                    self.types[i] = if perlin_idxs[j] == idx_cutoff {
                        CellType::Precinct(0)
                    } else {
                        if perlin_idxs[j] < idx_cutoff {
                            CellType::Precinct(-1)
                        } else {
                            CellType::Precinct(1)
                        }
                    };
                    j += 1;
                }
                _ => {}
            }
        }
    }

    // randomly assign sea or land to non-precinct cells
    fn randomize_others(&mut self, threshold: f64) {
        let perlin = perlin_noise::PerlinNoise::new();
        let num_cells = self.types.len();
        for i in 0..num_cells {
            match self.types[i] {
                CellType::Land | CellType::Sea => {
                    let p = perlin.get2d([
                        self.sites[i].x as f64 / self.width * 255.0 * 0.001 + 128.0,
                        self.sites[i].y as f64 / self.height * 255.0 * 0.001 + 128.0,
                    ]);
                    self.types[i] = if p > threshold {
                        CellType::Land
                    } else {
                        CellType::Sea
                    }
                }
                _ => {}
            }
        }
    }
}

fn create_map(num_precincts: usize) -> Result<GameMap, &'static str> {
    // manual parameters
    let w = 1280.0;
    let h = 720.0;
    let num_cells = num_precincts + 28;

    // create game map with all interior cells marked as precincts
    let mut game_map = GameMap::new(num_cells, w, h);

    // remove some random precincts near the edge
    let num_remove = game_map.num_precincts() as isize - num_precincts as isize - 1isize;
    for _ in 0..num_remove {
        game_map.remove_precinct_edge();
    }

    // remove precincts in small connected components
    game_map.remove_precinct_comp(num_precincts * 2 / 3);

    // make sure there are enough precincts left
    if game_map.num_precincts() < num_precincts {
        return Err("Map generation did not have enough precincts.");
    }

    // combine some neighboring precincts to bring the total down to the correct number of precincts
    while game_map.num_precincts() > num_precincts {
        game_map.combine_precincts()?;
    }

    // randomize precinct vote margins
    game_map.randomize_precincts();

    // randomize land/sea
    game_map.randomize_others(0.44);

    // return map and unused cells
    Ok(game_map)
}

fn spawn_polygon(
    commands: &mut Commands,
    vertices: &Vec<Vec2>,
    color: Color,
    outline_width: f32,
    z: f32,
) {
    let draw_mode = if outline_width > 0.0 {
        DrawMode::Outlined {
            fill_mode: FillMode::color(color),
            outline_mode: StrokeMode::new(Color::BLACK, outline_width),
        }
    } else {
        DrawMode::Fill(FillMode::color(color))
    };
    commands.spawn_bundle(GeometryBuilder::build_as(
        &shapes::Polygon {
            points: vertices.clone(),
            closed: true,
        },
        draw_mode,
        Transform::from_xyz(0.0, 0.0, z),
    ));
}

fn spawn_hills(commands: &mut Commands, vertices: &Vec<Vec2>, color: Color) {
    let mut rng = thread_rng();
    let mut center = Vec2::ZERO;
    for &v in vertices {
        center += v;
    }
    center /= vertices.len() as f32;
    let draw_mode = DrawMode::Fill(FillMode::color(color));
    let num_hills = rng.gen_range(0..3);
    let y_offset = rng.gen_range(-20.0..20.0);
    for _ in 0..num_hills {
        // hill vertices
        let w = rng.gen_range(50.0..80.0) as f32;
        let h = rng.gen_range(10.0..20.0) as f32;
        let v = vec![
            Vec2 { x: 0.5 * w, y: 0.0 },
            Vec2 { x: 0.0, y: h },
            Vec2 {
                x: -0.5 * w,
                y: 0.0,
            },
        ];
        // hill offset
        let x_offset = rng.gen_range(-20.0..20.0);
        commands.spawn_bundle(GeometryBuilder::build_as(
            &shapes::Polygon {
                points: v,
                closed: true,
            },
            draw_mode,
            Transform::from_xyz(center.x + x_offset, center.y + y_offset - 30.0, 2.0),
        ));
    }
}

fn spawn_waves(commands: &mut Commands, vertices: &Vec<Vec2>, color: Color) {
    let mut rng = thread_rng();
    let mut center = Vec2::ZERO;
    for &v in vertices {
        center += v;
    }
    center /= vertices.len() as f32;
    let mut path_builder = PathBuilder::new();
    path_builder.move_to(Vec2::ZERO);
    path_builder.line_to(Vec2 { x: 5.0, y: 0.0 });
    path_builder.line_to(Vec2 { x: 10.0, y: 5.0 });
    path_builder.line_to(Vec2 { x: 15.0, y: 8.0 });
    path_builder.line_to(Vec2 { x: 20.0, y: 8.0 });
    path_builder.line_to(Vec2 { x: 18.0, y: 4.0 });
    path_builder.line_to(Vec2 { x: 20.0, y: 2.0 });
    path_builder.line_to(Vec2 { x: 25.0, y: 0.0 });
    path_builder.line_to(Vec2 { x: 30.0, y: 0.0 });
    let line = path_builder.build();

    let x_offset = rng.gen_range(-20.0..20.0);
    let y_offset = rng.gen_range(-20.0..20.0);
    commands.spawn_bundle(GeometryBuilder::build_as(
        &line,
        DrawMode::Stroke(StrokeMode::new(color, 3.0)),
        Transform::from_xyz(center.x + x_offset - 15.0, center.y + y_offset, 3.0),
    ));
}

fn setup_system(
    mut commands: Commands,
    // asset_server: Res<AssetServer>,
    mut exit: EventWriter<bevy::app::AppExit>,
) {
    // manual settings
    let num_precincts = 25;

    let precinct_neutral_color = Color::rgb(0.6, 0.6, 0.6);
    let precinct_player_color = Color::rgb(0.784, 0.937, 0.682);
    let precinct_opponent_color = Color::rgb(0.375, 0.284, 0.371);
    let background_land_color = Color::rgb(0.58, 0.553, 0.498);
    let background_hill_color = Color::rgb(0.63, 0.603, 0.548);
    let background_sea_color = Color::rgb(0.522, 0.624, 0.659);
    let background_wave_color = Color::rgb(0.622, 0.724, 0.759);

    // create map
    let mut game_map_result: Result<GameMap, &'static str> = Err("Uninitialized game map");
    for _ in 0..5 {
        game_map_result = create_map(num_precincts);
        if let Err(e) = game_map_result {
            println!("Map generation failed: {:?}", e);
        } else {
            break;
        }
    }
    let game_map = match game_map_result {
        Ok(x) => x,
        Err(_) => {
            exit.send(bevy::app::AppExit);
            return;
        }
    };

    // camera
    commands.spawn_bundle(Camera2dBundle::default());

    // spawn non-precinct background polygons
    let background_land = game_map.land_vertices();
    let background_sea = game_map.sea_vertices();
    for vertices in background_land {
        spawn_polygon(&mut commands, &vertices, background_land_color, 0.0, 1.0);
        spawn_hills(&mut commands, &vertices, background_hill_color);
    }
    for vertices in background_sea {
        spawn_polygon(&mut commands, &vertices, background_sea_color, 0.0, 2.0);
        spawn_waves(&mut commands, &vertices, background_wave_color);
    }

    // spawn precinct polygons
    let num_precincts = game_map.num_precincts();
    let precinct_vertices = game_map.precinct_vertices();
    let precinct_margins = game_map.precinct_margins();
    // let precinct_sites = game_map.precinct_sites();
    // let precinct_indices = game_map.precinct_indices();
    // let font = asset_server.load("fonts/FiraSans-Bold.ttf");
    // let text_style = TextStyle {
    //     font,
    //     font_size: 60.0,
    //     color: Color::BLACK,
    // };
    for i in 0..num_precincts {
        let fill_color = match precinct_margins[i] {
            1 => precinct_player_color,
            -1 => precinct_opponent_color,
            _ => precinct_neutral_color,
        };
        spawn_polygon(&mut commands, &precinct_vertices[i], fill_color, 2.0, 4.0);
        // commands.spawn_bundle(Text2dBundle {
        //     text: Text::from_section(precinct_indices[i].to_string(), text_style.clone())
        //         .with_alignment(TextAlignment::CENTER),
        //     transform: Transform::from_xyz(
        //         precinct_sites[i].x, // - 60.0 / 2.0,
        //         precinct_sites[i].y, // + 60.0 / 2.0,
        //         1.0,
        //     ),
        //     ..default()
        // });
    }
}

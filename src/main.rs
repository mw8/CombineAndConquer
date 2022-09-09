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
        .add_system(update_system)
        .run();
}

#[derive(Clone, Copy)]
enum CellType {
    Empty,         // ignore this cell
    Precinct(i32), // precinct with given vote margin
    Land,          // non-precinct land cell
    Sea,           // non-precinct sea cell
}

struct GameMap {
    width: f64,
    height: f64,
    types: Vec<CellType>,       // type of each cell
    sites: Vec<Vec2>,           // point near the center of each cell
    dist2: Vec<f32>, // sq. radius of bounding circle (centered at each site) of each cell
    vertices: Vec<Vec<Vec2>>, // vertices of each cell
    neighbors: Vec<Vec<usize>>, // neighbors of each cell
    districts: Vec<Vec<usize>>, // collections of cells to form districts
}

impl GameMap {
    // create a new random game map
    fn new(num_cells: usize, width: f64, height: f64) -> GameMap {
        // random starting sites for cells
        let mut rng = thread_rng();
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
            .set_lloyd_relaxation_iterations(4)
            .build()
            .unwrap();
        // extract voronoi cell properties
        let mut types: Vec<CellType> = Vec::with_capacity(num_cells);
        let mut sites: Vec<Vec2> = Vec::with_capacity(num_cells);
        let mut dist2: Vec<f32> = Vec::with_capacity(num_cells);
        let mut vertices: Vec<Vec<Vec2>> = Vec::with_capacity(num_cells);
        for cell in voronoi.iter_cells() {
            let site = cell.site_position();
            let cell_site = Vec2::new(site.x as f32, site.y as f32);
            sites.push(cell_site);
            let mut cell_vertices: Vec<Vec2> = Vec::with_capacity(cell.iter_vertices().count());
            let mut max_dist2 = 0.0f32;
            for vertex in cell.iter_vertices() {
                let v = Vec2::new(vertex.x as f32, vertex.y as f32);
                let d2 = cell_site.distance_squared(v);
                if max_dist2 < d2 {
                    max_dist2 = d2;
                }
                cell_vertices.push(v);
            }
            dist2.push(max_dist2 + 0.1);
            vertices.push(cell_vertices);
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
            width: width,
            height: height,
            types: types,
            sites: sites,
            dist2: dist2,
            vertices: vertices,
            neighbors: neighbors,
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
    fn remove_precinct_edge(&mut self, rng: &mut ThreadRng) {
        let num_cells = self.num_cells();
        let mut i0 = rng.gen_range(0..num_cells);
        let mut j0 = 0;
        while j0 < num_cells {
            if matches!(self.types[i0], CellType::Land) {
                let neighbors = &self.neighbors[i0];
                let num_neighbors = neighbors.len();
                let mut i1 = rng.gen_range(0..num_neighbors);
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
        let num_cells = self.num_cells();
        // remove i0 from i1 neighbors
        let mut nb = self.neighbors[i1].clone();
        if let Some(i0_pos) = nb.iter().position(|&x| x == i0) {
            nb.remove(i0_pos);
        } else {
            return Err("Cannot merge cells because they are not neighbors");
        }
        // remove i1 from i0 neighbors
        if let Some(i1_pos) = self.neighbors[i0].iter().position(|&x| x == i1) {
            self.neighbors[i0].remove(i1_pos);
        } else {
            return Err("Cannot merge cells because they are not neighbors");
        }
        // merge neighbors
        self.neighbors[i0].append(&mut nb);
        self.neighbors[i1].clear();
        // find all instances of i1 in other neighbor lists and change them to i0
        for i in 0..num_cells {
            if i != i0 && i != i1 {
                for ni in self.neighbors[i].iter_mut() {
                    if *ni == i1 {
                        *ni = i0;
                    }
                }
            }
        }
        // merge vertices and normals
        let merged_vertices =
            merge_cell_vertices(self.vertices[i0].clone(), self.vertices[i1].clone())?;
        self.vertices[i0] = merged_vertices;
        self.vertices[i1].clear();
        // update types and sites
        self.types[i1] = CellType::Empty;
        self.sites[i0] = (self.sites[i0] + self.sites[i1]) * 0.5;
        Ok(())
    }

    // combine two random neighboring precincts
    fn combine_precincts(&mut self, rng: &mut ThreadRng) -> Result<(), &'static str> {
        // find random precinct
        let num_cells = self.num_cells();
        let mut i0 = rng.gen_range(0..num_cells);
        let mut j = 0;
        while j < num_cells && !matches!(self.types[i0], CellType::Precinct(_)) {
            i0 = (i0 + 1) % num_cells;
            j += 1;
        }
        if j == num_cells {
            return Err("Unable to find any precincts to combine.");
        }
        // find a neighboring precinct
        let neighbors = &self.neighbors[i0];
        let num_neighbors = neighbors.len();
        let mut i1 = rng.gen_range(0..num_neighbors);
        j = 0;
        while j < num_neighbors && !matches!(self.types[neighbors[i1]], CellType::Precinct(_)) {
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

    // test if a given cell (index i) contains a given point
    fn within_cell(&self, i: usize, p: Vec2) -> bool {
        assert!(i < self.num_cells());
        let n = self.vertices[i].len();
        for j in 0..n {
            let v = self.vertices[i][j];
            let d = self.vertices[i][(j + 1) % n] - v;
            if d.y * (p.x - v.x) - d.x * (p.y - v.y) < 0.0 {
                return false;
            }
        }
        true
    }

    // select a cell given some point in the map
    fn select_cell(&self, p: Vec2) -> Result<usize, &'static str> {
        if p.x < -0.5 * self.width as f32
            || p.x > 0.5 * self.width as f32
            || p.y < -0.5 * self.height as f32
            || p.y > 0.5 * self.height as f32
        {
            return Err("Attempted to select cell with coordinates outside the game map.");
        }
        let num_cells = self.num_cells();
        for i in 0..num_cells {
            if p.distance_squared(self.sites[i]) < self.dist2[i] && self.within_cell(i, p) {
                return Ok(i);
            }
        }
        Err("Could not find any cell containing those coordinates")
    }
}

fn setup_system(mut commands: Commands, mut exit: EventWriter<bevy::app::AppExit>) {
    // camera
    commands.spawn_bundle(Camera2dBundle::default());
    // map
    let game_map = match spawn_game_map(&mut commands, 25) {
        Ok(x) => x,
        Err(e) => {
            println!("{:?}", e);
            exit.send(bevy::app::AppExit);
            return;
        }
    };
    commands.insert_resource(game_map);
}

fn update_system(
    mut commands: Commands,
    buttons: Res<Input<MouseButton>>,
    windows: Res<Windows>,
    mut game_map: ResMut<GameMap>,
    district_query: Query<Entity, With<DistrictTag>>,
) {
    let win_size = Vec2::new(1280.0, 720.0);
    if buttons.just_pressed(MouseButton::Left) {
        let window = windows.get_primary().unwrap();
        if let Some(pos_win) = window.cursor_position() {
            // compute click position
            let pos = pos_win - win_size * 0.5;
            println!("Mouse button pressed at location: {:?}!", pos);
            // find cell selected
            match game_map.select_cell(pos) {
                Ok(idx) => {
                    if game_map.districts.is_empty() {
                        game_map.districts.push(Vec::new());
                    }
                    game_map.districts[0].push(idx);
                }
                Err(e) => {
                    println!("Error selecting cell: {:?}", e);
                }
            }
            // despawn all current districts
            for district_entity in &district_query {
                commands.entity(district_entity).despawn();
            }
            // spawn updated districts
            spawn_districts(&mut commands, &game_map, Color::ORANGE);
        }
    }
}

fn create_map(num_precincts: usize) -> Result<GameMap, &'static str> {
    // manual parameters
    let w = 1280.0;
    let h = 720.0;

    let num_cells = num_precincts + 30;
    let mut rng = thread_rng();

    // create game map with all interior cells marked as precincts
    let mut game_map = GameMap::new(num_cells, w, h);

    // remove some random precincts near the edge
    let num_remove = game_map.num_precincts() as isize - num_precincts as isize - 2isize;
    for _ in 0..num_remove {
        game_map.remove_precinct_edge(&mut rng);
    }

    // remove precincts in small connected components
    game_map.remove_precinct_comp(num_precincts * 2 / 3);

    // make sure there are enough precincts left
    if game_map.num_precincts() < num_precincts {
        return Err("Map generation did not have enough precincts.");
    }

    // combine some neighboring precincts to bring the total down to the correct number of precincts
    while game_map.num_precincts() > num_precincts {
        game_map.combine_precincts(&mut rng)?;
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
    line_width: f32,
    x_offset: f32,
    y_offset: f32,
    z_offset: f32,
) {
    let draw_mode = if line_width > 0.0 {
        DrawMode::Outlined {
            fill_mode: FillMode::color(color),
            outline_mode: StrokeMode::new(Color::BLACK, line_width),
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
        Transform::from_xyz(x_offset, y_offset, z_offset),
    ));
}

fn spawn_path(
    commands: &mut Commands,
    vertices: &Vec<Vec2>,
    color: Color,
    line_width: f32,
    closed: bool,
    x_offset: f32,
    y_offset: f32,
    z_offset: f32,
) {
    let n = vertices.len();
    if n < 2 {
        return;
    }
    let mut path_builder = PathBuilder::new();
    path_builder.move_to(vertices[0]);
    for i in 1..n {
        path_builder.line_to(vertices[i]);
    }
    if closed {
        path_builder.line_to(vertices[0]);
    }
    let path = path_builder.build();
    commands.spawn_bundle(GeometryBuilder::build_as(
        &path,
        DrawMode::Stroke(StrokeMode::new(color, line_width)),
        Transform::from_xyz(x_offset, y_offset, z_offset),
    ));
}

fn spawn_hills(commands: &mut Commands, cell_vertices: &Vec<Vec2>, color: Color) {
    let mut rng = thread_rng();
    let mut center = Vec2::ZERO;
    for &v in cell_vertices {
        center += v;
    }
    center /= cell_vertices.len() as f32;
    let num_hills = rng.gen_range(0..3);
    let y_offset = rng.gen_range(-20.0..20.0) - 30.0;
    for _ in 0..num_hills {
        let x_offset = rng.gen_range(-20.0..20.0);
        let w = rng.gen_range(50.0..80.0) as f32;
        let h = rng.gen_range(10.0..20.0) as f32;
        let hill_vertices = vec![
            Vec2 { x: 0.5 * w, y: 0.0 },
            Vec2 { x: 0.0, y: h },
            Vec2 {
                x: -0.5 * w,
                y: 0.0,
            },
        ];
        spawn_polygon(
            commands,
            &hill_vertices,
            color,
            0.0,
            center.x + x_offset,
            center.y + y_offset,
            2.0,
        )
    }
}

fn spawn_waves(commands: &mut Commands, cell_vertices: &Vec<Vec2>, color: Color) {
    let mut rng = thread_rng();
    let mut center = Vec2::ZERO;
    for &v in cell_vertices {
        center += v;
    }
    center /= cell_vertices.len() as f32;
    let wave_vertices = vec![
        Vec2::ZERO,
        Vec2 { x: 5.0, y: 0.0 },
        Vec2 { x: 10.0, y: 5.0 },
        Vec2 { x: 15.0, y: 8.0 },
        Vec2 { x: 20.0, y: 8.0 },
        Vec2 { x: 18.0, y: 4.0 },
        Vec2 { x: 20.0, y: 2.0 },
        Vec2 { x: 25.0, y: 0.0 },
        Vec2 { x: 30.0, y: 0.0 },
    ];
    let x_offset = rng.gen_range(-20.0..20.0);
    let y_offset = rng.gen_range(-20.0..20.0);
    spawn_path(
        commands,
        &wave_vertices,
        color,
        3.0,
        false,
        center.x + x_offset - 15.0,
        center.y + y_offset,
        3.0,
    );
}

#[derive(Component)]
struct DistrictTag;

fn spawn_districts(commands: &mut Commands, game_map: &GameMap, color: Color) {
    let line_width = 12.0;
    let line_offset = 8.0;
    let draw_mode = DrawMode::Stroke(StrokeMode::new(color, line_width));
    for district in &game_map.districts {
        let district_size = district.len();
        if district_size > 0 {
            // make district vertices
            let mut v0 = game_map.vertices[district[0]].clone();
            for i in 1..district_size {
                let v1 = game_map.vertices[district[i]].clone();
                match merge_cell_vertices(v0.clone(), v1) {
                    Ok(v2) => {
                        v0 = v2;
                    }
                    Err(e) => {
                        println!("District generation failed: {:?}", e);
                        return;
                    }
                }
            }
            // go through all vertices and move them to the intersection point of
            // the offset lines going parallel to the vertices of the district
            match offset_cell_vertices(&v0, line_offset) {
                Ok(v) => {
                    // spawn district polygon
                    commands
                        .spawn_bundle(GeometryBuilder::build_as(
                            &shapes::Polygon {
                                points: v,
                                closed: true,
                            },
                            draw_mode,
                            Transform::from_xyz(0.0, 0.0, 5.0),
                        ))
                        .insert(DistrictTag);
                }
                Err(e) => {
                    println!("District generation failed: {:?}", e);
                    return;
                }
            }
        }
    }
}

fn spawn_game_map(commands: &mut Commands, num_precincts: usize) -> Result<GameMap, &'static str> {
    // colors
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
    let game_map = game_map_result?;
    // spawn non-precinct background polygons
    let background_land = game_map.land_vertices();
    let background_sea = game_map.sea_vertices();
    for vertices in background_land {
        spawn_polygon(
            commands,
            &vertices,
            background_land_color,
            0.0,
            0.0,
            0.0,
            1.0,
        );
        spawn_hills(commands, &vertices, background_hill_color);
    }
    for vertices in background_sea {
        spawn_polygon(
            commands,
            &vertices,
            background_sea_color,
            0.0,
            0.0,
            0.0,
            2.0,
        );
        spawn_waves(commands, &vertices, background_wave_color);
    }
    // spawn precinct polygons
    let num_precincts = game_map.num_precincts();
    let precinct_vertices = game_map.precinct_vertices();
    let precinct_margins = game_map.precinct_margins();
    for i in 0..num_precincts {
        let fill_color = match precinct_margins[i] {
            1 => precinct_player_color,
            -1 => precinct_opponent_color,
            _ => precinct_neutral_color,
        };
        spawn_polygon(
            commands,
            &precinct_vertices[i],
            fill_color,
            3.0,
            0.0,
            0.0,
            4.0,
        );
    }
    // return game map
    Ok(game_map)
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

// circularly shift a vector so that the element at a given vertex moves to the beginning
fn circ_shift_vec<T>(v: &mut Vec<T>, i: usize) {
    let mut v0 = v.split_off(i);
    v0.append(v);
    *v = v0;
}

// merge vertices for two neighboring cells
fn merge_cell_vertices(mut v0: Vec<Vec2>, mut v1: Vec<Vec2>) -> Result<Vec<Vec2>, &'static str> {
    // find a vertex u in common between the two cells
    let mut u = Vec2 {
        x: std::f32::NAN,
        y: std::f32::NAN,
    };
    let mut found = false;
    for &u0 in &v0 {
        for &u1 in &v1 {
            if u0 == u1 {
                u = u0;
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
    let mut u0 = v0.iter().position(|&x| x == u).unwrap();
    let mut u1 = v1.iter().position(|&x| x == u).unwrap();
    assert!(v0[u0] == v1[u1]);
    let mut w0 = u0;
    let mut w1 = u1;
    // move u as far forward/backward as possible in v0/v1
    let mut j0 = (u0 + v0.len() - 1) % v0.len();
    let mut j1 = (u1 + v1.len() + 1) % v1.len();
    while v0[j0] == v1[j1] {
        u0 = j0;
        u1 = j1;
        j0 = (j0 + v0.len() - 1) % v0.len();
        j1 = (j1 + v1.len() + 1) % v1.len();
    }
    assert!(v0[u0] == v1[u1]);
    u = v0[u0];
    // move w as far backward/forward as possible in v0/v1
    j0 = (w0 + v0.len() + 1) % v0.len();
    j1 = (w1 + v1.len() - 1) % v1.len();
    while v0[j0] == v1[j1] {
        w0 = j0;
        w1 = j1;
        j0 = (j0 + v0.len() + 1) % v0.len();
        j1 = (j1 + v1.len() - 1) % v1.len();
    }
    assert!(v0[w0] == v1[w1]);
    let w = v0[w0];
    if u == w {
        return Err("Cannot merge cells because they have only 1 vertex in common");
    }
    // circularly shift v0/n0 so that v is index 0 and v1/n1 so that w is index 0
    circ_shift_vec(&mut v0, u0);
    assert!(v0[0] == u);
    circ_shift_vec(&mut v1, w1);
    assert!(v1[0] == w);
    w0 = ((w0 + v0.len()) - u0) % v0.len();
    u1 = ((u1 + v1.len()) - w1) % v1.len();
    assert!(v0[w0] == w);
    assert!(v1[u1] == u);
    // merge cell vertices
    let mut v01 = v0.split_off(w0);
    let mut v11 = v1.split_off(u1);
    assert!(v0.len() > 0);
    assert!(v1.len() > 0);
    assert!(v01.len() > 0);
    assert!(v11.len() > 0);
    v01.append(&mut v11);
    Ok(v01)
}

fn offset_cell_vertices(vertices: &Vec<Vec2>, offset: f32) -> Result<Vec<Vec2>, &'static str> {
    let n = vertices.len();
    if n < 3 {
        return Err("Need 3 or more vertices to offset cell.");
    }
    let mut w = Vec::new();
    for i in 0..n {
        // get current vertex v1 and neighboring vertices v0 and v2
        let v0 = vertices[(i + n - 1) % n];
        let v1 = vertices[i];
        let v2 = vertices[(i + 1) % n];
        // compute difference vectors between neighboring vertices
        let d0 = v1 - v0;
        let d1 = v2 - v1;
        let p0 = Vec2::new(d0.y, -d0.x).normalize();
        let p1 = Vec2::new(d1.y, -d1.x).normalize();
        // solve for
        let det = d1.x * d0.y - d0.x * d1.y;
        if (d0.x * d1.y - d1.x * d0.y).abs() < 1e-3 || det.abs() < 1e-3 {
            continue;
        }
        let b0 = d0.x + (p1.x - p0.x) * offset;
        let b1 = d0.y + (p1.y - p0.y) * offset;
        let s = (d1.x * b1 - d1.y * b0) / det;
        w.push(Vec2::new(
            v0.x + p0.x * offset + d0.x * s,
            v0.y + p0.y * offset + d0.y * s,
        ));
    }
    Ok(w)
}

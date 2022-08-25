use bevy::prelude::*;
use bevy_prototype_lyon::prelude::{tess::path::Polygon, *};
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

fn setup_system(mut commands: Commands) {
    let n = 75;
    let w = 1250.0;
    let h = 690.0;
    // voronoi sites
    let mut rng = thread_rng();
    let mut points0 = Vec::<Point>::with_capacity(n);
    for _ in 0..n {
        let x = rng.gen_range(-0.5 * w..0.5 * w);
        let y = rng.gen_range(-0.5 * h..0.5 * h);
        points0.push(Point { x: x, y: y })
    }
    println!("{:?}", points0);
    // builds a voronoi diagram from the set of sites above, bounded by a square of size 4
    let voronoi = VoronoiBuilder::default()
        .set_sites(points0)
        .set_bounding_box(BoundingBox::new_centered(w, h))
        .set_lloyd_relaxation_iterations(1)
        .build()
        .unwrap();
    // display polygons
    let x0 = (0.5 * w - 1e-3) as f32;
    let y0 = (0.5 * h - 1e-3) as f32;
    let mut points1: Vec<Vec2> = Vec::new();
    commands.spawn_bundle(Camera2dBundle::default());
    for cell in voronoi.iter_cells() {
        println!("{:?}", cell);
        // extract points
        points1.clear();
        for vertex in cell.iter_vertices() {
            points1.push(Vec2::new(vertex.x as f32, vertex.y as f32));
        }
        // determine if boundary cell
        let mut boundary = false;
        for p in points1.iter() {
            println!("{:?}", p);
            if p.x < -x0 || p.x > x0 || p.y < -y0 || p.y > y0 {
                boundary = true;
                break;
            }
        }
        // create polygon with points
        let shape = shapes::Polygon {
            points: points1.clone(),
            closed: true,
        };
        // spawn polygon
        let fill_mode = if boundary {
            FillMode::color(Color::CYAN)
        } else {
            FillMode::color(Color::GRAY)
        };
        commands.spawn_bundle(GeometryBuilder::build_as(
            &shape,
            DrawMode::Outlined {
                fill_mode: fill_mode,
                outline_mode: StrokeMode::new(Color::BLACK, 5.0),
            },
            Transform::default(),
        ));
    }
}

//     // inspect cells through iterators
//     my_voronoi.iter_cells().for_each(|cell| {
//         println!(
//             "Vertices of cell: {:?}",
//             cell.iter_vertices().collect::<Vec<&Point>>()
//         )
//     });

//     // or probe cells individually
//     let my_cell = my_voronoi.cell(1);
//     println!(
//         "Second cell has site {:?}, voronoi vertices {:?} and delaunay triangles {:?}",
//         my_cell.site_position(),
//         my_cell.iter_vertices().collect::<Vec<&Point>>(),
//         my_cell.triangles().iter().collect::<Vec<&usize>>()
//     );
// }

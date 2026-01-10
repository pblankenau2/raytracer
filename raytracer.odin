package raytracer

import "core:sort"
import "core:fmt"
import "core:os"
import "core:math/linalg"
import "core:math"
import "core:strings"
import "core:slice"
import "core:testing"

make_vec3 :: proc(x:f64,y:f64,z:f64) -> [4]f64 {
	return [4]f64{x,y,z,0}
}


make_pnt3 :: proc(x:f64,y:f64,z:f64) -> [4]f64 {
	return [4]f64{x,y,z,1}
}

Color :: [3]f64

Canvas :: struct {
	width: int,
	height: int,
	pixels: [dynamic]Color
}

Ray :: struct {
	origin: [4]f64, // Point
	direction: [4]f64 // Vec
}

Sphere :: struct {
	obj_id: i32,
	transform: matrix[4,4]f64,
	material: Material
}

make_sphere :: proc(
	obj_id: i32,
	transform: matrix[4,4]f64 = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1},
	material : Material = DefaultMaterial
) -> Sphere {
	return Sphere{obj_id=obj_id, transform=transform, material=material}
}

Intersection :: struct {
	t: f64,
	object: Sphere, // TODO: make it a generic Object type.
}

LightPoint :: struct {
	position: [4]f64,
	intensity: Color
}

Material :: struct {
	color: Color,
	ambient: f64,
	diffuse: f64,
	specular: f64,
	shininess: f64
}

DefaultMaterial :: Material{Color{1,1,1}, 0.1, 0.9, 0.9, 200.0}

World :: struct {
	// We might want to have struct of arrays here.  Array of spheres, arrays other objects
	light: LightPoint,
	spheres: []Sphere
}


DefaultWorld := World{
	LightPoint{make_pnt3(-10,10,-10),Color{1,1,1}},
	{
		make_sphere(0, material=Material{color=Color{0.8,1.0,0.6}, ambient=0.1, diffuse=0.7, specular=0.2, shininess=200}),
		make_sphere(1, transform=linalg.matrix4_scale([3]f64{0.5,0.5,0.5}))
	}
}

intersect_world :: proc(world: World, ray: Ray) -> []Intersection {
	intersections : [dynamic]Intersection // heap allocated
	for s, i in world.spheres {
		i1, i2, okay := intersect(ray, s)
		if okay {
			append(&intersections, i1, i2)
		}
	}
	slice.sort_by(intersections[:], proc(i, j: Intersection) -> bool {
        return i.t < j.t
    })

	return intersections[:] // need the [:] to make it a slice?
}

@(test)
test_intersect_world :: proc(t: ^testing.T) {

	r := Ray{make_pnt3(0,0,-5), make_vec3(0,0,1)}
	xs := intersect_world(DefaultWorld, r)

	testing.expect(t, len(xs) == 4)
	testing.expect(t, xs[0].t == 4)
	testing.expect(t, xs[1].t == 4.5)
	testing.expect(t, xs[2].t == 5.5)
	testing.expect(t, xs[3].t == 6)
}

PreComputations :: struct {
	t: f64,
	object: Sphere,
	point: [4]f64,
	eyev: [4]f64,
	normalv: [4]f64,
	inside: bool
}

prepare_computations :: proc(intersection: Intersection, ray: Ray) -> PreComputations {
	point := position(ray, intersection.t)
	normalv := normal_at(intersection.object, point)
	eyev := -ray.direction

	inside := false
	if linalg.dot(normalv, eyev) < 0 {
		inside = true
		normalv = -normalv
	}

	return PreComputations{
		t=intersection.t,
		object=intersection.object,
		point=point,
		eyev=eyev,
		normalv=normalv,
		inside=inside
	}
}

shade_hit :: proc(world: World, comps: PreComputations) -> Color {
	return lighting(comps.object.material, world.light, comps.point, comps.eyev, comps.normalv)
}

@(test)
test_shade_hit :: proc(t: ^testing.T) {
	r := Ray{make_pnt3(0,0,-5), make_vec3(0,0,1)}
	i := Intersection{4, DefaultWorld.spheres[0]}
	c := shade_hit(DefaultWorld, prepare_computations(i, r))
	testing.expect(t, linalg.vector_length(c - Color{0.38066, 0.47583, 0.2855}) < f64(0.001))


	w := DefaultWorld // TODO: Is this a copy? Try in debugger to find out...
	w.light = LightPoint{make_pnt3(0, 0.25, 0), Color{1,1,1}}
	r = Ray{make_pnt3(0, 0, 0), make_vec3(0, 0, 1)}
	i = Intersection{0.5, w.spheres[1]}

	c = shade_hit(w, prepare_computations(i, r))
	testing.expect(t, linalg.vector_length(c - Color{0.90498, 0.90498, 0.90498}) < f64(0.001))
}

lighting :: proc(material: Material, light: LightPoint, point: [4]f64, eyev: [4]f64, normalv: [4]f64) -> Color {
	effective_color := material.color * light.intensity
	lightv := linalg.normalize(light.position - point)
	ambient := effective_color * material.ambient
	light_dot_normal := linalg.dot(lightv, normalv)

	diffuse := Color{0,0,0}
	specular := Color{0,0,0}
	if light_dot_normal < 0 {
		diffuse = Color{0,0,0}
		specular = Color{0,0,0}
	} else {
		diffuse = effective_color * material.diffuse * light_dot_normal
		reflectv := reflect(-lightv, normalv)
		reflect_dot_eye := linalg.dot(reflectv, eyev)
		if reflect_dot_eye > 0 {
			factor := math.pow(reflect_dot_eye, material.shininess)
			specular = light.intensity * material.specular * factor
		}
	}
	// fmt.println(ambient, diffuse, specular)
	return ambient + diffuse + specular
}


// TODO: return [2]f64 slice?
intersect :: proc(ray: Ray, object: Sphere) -> (Intersection, Intersection, bool) { 

	// Need to transform the ray before calculating the intersection.
	new_ray := transform(ray, linalg.inverse(object.transform))
	sphere_to_ray := new_ray.origin - make_pnt3(0.0,0.0,0.0)
	a := linalg.dot(new_ray.direction, new_ray.direction)
	b := 2 * linalg.dot(new_ray.direction, sphere_to_ray)
	c := linalg.dot(sphere_to_ray, sphere_to_ray) - 1

	discriminant := math.pow(b,2) - 4*a*c
	if (discriminant < 0.0) {return Intersection{0.0, object}, Intersection{0.0, object}, false} // TODO: return nil bad?

	t1 := (-b - math.sqrt(discriminant)) / (2 * a)
	t2 := (-b + math.sqrt(discriminant)) / (2 * a)

	// if ray is tangent to sphere then return the same intersection twice
	return Intersection{t1, object}, Intersection{t2, object}, true
}

// TODO: is this function needed?  Can we combine with intersection?
hit :: proc(intersections: []Intersection) -> (Intersection, bool) {
	
	smallest := intersections[0]
	for i in intersections {
		// Remove negative t intersections
		if (i.t < 0) {continue}

		// Hit is smallest non negative intersection
		if (i.t < smallest.t) {smallest = i}
	}
	if (smallest.t < 0) {
		return smallest, false
	}
	return smallest, true
}

position :: proc(ray: Ray, t: f64) -> [4]f64 {
	return ray.origin + ray.direction * t
}

set_color :: proc(canvas: ^Canvas, row: int, col: int, color: Color) {
	clipped_color := Color{min(color.r, 1), min(color.g, 1), min(color.b, 1)}
	if (row < canvas.height) & (col < canvas.width) & (row >= 0) & (col >= 0) {
		canvas.pixels[canvas.width*row + col] = clipped_color
	}

}

make_canvas :: proc(width: int, height: int) -> Canvas {
	pix: [dynamic]Color

	for i in 0..<(width*height) {
		append(&pix, Color{0,0,0})
	}
	// return Canvas{width, height, pix}
	return Canvas{width, height, pix}
}

canvas_to_ppm :: proc(canvas: Canvas, filename: string) {
	handle, err := os.open("test.ppm", mode=(os.O_CREATE|os.O_TRUNC))
	if err!=nil {
		fmt.println("error!")
	}

	// os.write_string(handle, "P3\n")
	// os.write_string(handle, fmt.tprintf("%v %v\n", canvas.width, canvas.height))
	// os.write_string(handle, "255\n")

	builder := strings.builder_make(context.temp_allocator)
	strings.write_string(&builder, "P3\n")
	strings.write_string(&builder, fmt.tprintf("%v %v\n", canvas.width, canvas.height))
	strings.write_string(&builder, "255\n")

	counter := 0
	str_template := "%v %v %v "
	str_len := len(str_template)
	for v in canvas.pixels {
		r, g, b := int(v.r * 255), int(v.g * 255), int(v.b * 255)
		// r = int(v.r * 255)
		if (counter >= 70) {
			strings.write_string(&builder, "\n")
			strings.write_string(&builder, fmt.tprintf(str_template, r, g, b))
			counter = str_len
		} else {
			strings.write_string(&builder, fmt.tprintf(str_template, r, g, b))
			counter = counter + str_len
		}
	}
	str := strings.to_string(builder)
	os.write_string(handle, str)
}

transform :: proc(r: Ray, m: matrix[4,4]f64) -> Ray {
	return Ray{m * r.origin, m * r.direction}
}

normal_at :: proc(s: Sphere, p: [4]f64) -> [4]f64 {
	object_point := linalg.inverse(s.transform) * p
	object_normal := object_point - [4]f64{0.0,0.0,0.0,0.0}
	world_normal := linalg.transpose(linalg.inverse(s.transform)) * object_normal
	world_normal.w = 0.0
	return linalg.normalize(world_normal)
}

reflect :: proc(in_: [4]f64, normal: [4]f64) -> [4]f64 {
	return in_ - normal * 2 * linalg.dot(in_, normal)
}

main :: proc() {
	origin := make_pnt3(0,0,-5)
	s := make_sphere(0)
	// s.transform = linalg.matrix4_scale([3]f64{1,0.5,1}) // TODO: No need for set_transform function.
	s.material.color = Color{1, 0.2, 1}
	s.material.shininess = 100

	light_position := make_pnt3(-10, 10, -10)
	light_color := Color{1, 1, 1}
	light := LightPoint{light_position, light_color} 


	// Cast one ray per pixel?
	canv_rows := 500
	canv_cols := 500
	canv_world_size := 5.0
	canv_z := 5.0
	pixel_size := canv_world_size / f64(canv_rows)
	top_y := canv_world_size / 2
	left_x := -canv_world_size / 2
	canv := make_canvas(canv_cols, canv_rows)

	for row in 0..<canv_rows {
		world_y := top_y - f64(row)*pixel_size
		for col in 0..<canv_cols {
			// Canvas passes through sphere origin, pixels 100 pixels is 1 unit world space
			// So the top left pixel in world space is -25, -25, 0? the next pixel to the right is -25, -24.9, 0
			world_x := left_x + f64(col)*pixel_size

			pos := make_vec3(world_x,world_y,canv_z)
			ray := Ray{origin, linalg.normalize(pos - origin)}
			
			i1, i2, okay := intersect(ray, s)
			hit : Intersection
			if (i1.t < i2.t) & (i1.t > 0) {
				hit = i1
			} else if (i2.t < i1.t) & (i2.t > 0) {
				hit = i2
			} else {
				okay = false
			}

			if okay {
				point := position(ray, hit.t)
				normal := normal_at(hit.object, point)
				eye := -ray.direction
				color := lighting(hit.object.material, light, point, eye, normal)
				set_color(&canv, row, col, color)
			} else {
				set_color(&canv, row, col, Color{0.0, 0.0, 0.0})
			}
		}
	}

	canvas_to_ppm(canv, "test.ppm")

	r := Ray{make_pnt3(0,0,-5), make_vec3(0,0,1)}
	i := Intersection{4, DefaultWorld.spheres[0]}
	// fmt.println(shade_hit(DefaultWorld, prepare_computations(i, r)))
	comps := prepare_computations(i, r)
	c := shade_hit(DefaultWorld, comps)

	fmt.println(c)
	// testing.expect(t, shade_hit(DefaultWorld, prepare_computations(i, r)) == Color{0.38066, 0.47583, 0.2855})


	// w := DefaultWorld // TODO: Is this a copy? Try in debugger to find out...
	// w.light = LightPoint{make_pnt3(0, 0.25, 0), Color{1,1,1}}
	// r = Ray{make_pnt3(0, 0, 0), make_vec3(0, 0, 1)}
	// i = Intersection{0.5, w.spheres[1]}
	// fmt.println(shade_hit(w, prepare_computations(i, r)))
	// testing.expect(t, shade_hit(w, prepare_computations(i, r)) == Color{0.90498, 0.90498, 0.90498})
}

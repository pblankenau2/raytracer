package raytracer

import "core:sort"
import "core:fmt"
import "core:os"
import "core:math/linalg"
import "core:math"
import "core:strings"
import "core:slice"
import "core:testing"
import "vendor:stb/image"

EPSILON :: 0.0001

ColorU8 :: [4]u8

canvas_to_bmp :: proc(canv: Canvas) -> []ColorU8 {
	pixels : [dynamic]ColorU8
	for p in canv.pixels {
		append(&pixels, ColorU8{u8(p.r*255), u8(p.g*255), u8(p.b*255), 255})
	}
	return pixels[:]
}

make_vec3 :: proc(x:f64,y:f64,z:f64) -> [4]f64 {
	return [4]f64{x,y,z,0}
}


make_pnt3 :: proc(x:f64,y:f64,z:f64) -> [4]f64 {
	return [4]f64{x,y,z,1}
}

Color :: [3]f64

Canvas :: struct {
	width: i32,
	height: i32,
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

// NOTE: the canvas is always 1 world unit from the camera.
Camera :: struct {
	hsize: i32,
	vsize: i32,
	field_of_view: f64, //= (math.PI / 2.0)
	transform: matrix[4,4]f64,
	pixel_size: f64,
	half_width: f64,
	half_height: f64
}

make_camera :: proc(
	hsize: i32 = 200,
	vsize: i32 = 200,
	field_of_view: f64 = (math.PI / 2.0),
	transform: matrix[4,4]f64 = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}
) -> Camera { // TODO: would you ever want to make only have the Camera struct defined within a proc?
	half_view := math.tan(field_of_view / 2)
	aspect := f64(hsize) / f64(vsize)

	half_width : f64
	half_height : f64
	if aspect >= 1 {
		half_width = half_view
		half_height = half_view / aspect
	} else {
		half_width = half_view * aspect
		half_height = half_view
	}
	pixel_size := (half_width * 2) / f64(hsize)

	return Camera{
		hsize,
		vsize,
		field_of_view,
		transform,
		pixel_size,
		half_width,
		half_height
	}
}

ray_for_pixel :: proc(camera: Camera, px: i32, py: i32) -> Ray {
	// the offset from the edge of the canvas to the pixel's center
	xoffset := (f64(px) + 0.5) * camera.pixel_size
	yoffset := (f64(py) + 0.5) * camera.pixel_size

	// The untransformed coordinates of the pixel in world space.
	// Remember that the camera looks toward -z, so +x is to the left.
	world_x := camera.half_width - xoffset
	world_y := camera.half_height - yoffset

	// using the camera matrix, transform the canvas point and the origin,
	// and then compute the ray's direction vector. (remember that the canvas is at z=-1)
	pixel := linalg.inverse(camera.transform) * make_pnt3(world_x, world_y, -1)
	origin := linalg.inverse(camera.transform) * make_pnt3(0,0,0)
	direction := linalg.normalize(pixel - origin)
	return Ray{origin, direction}
}

render :: proc(camera: Camera, world: World) -> Canvas {
	canv := Canvas{
		width=camera.hsize,
		height=camera.vsize,
		pixels=make([dynamic]Color, camera.hsize*camera.vsize) // TODO: preallocated dynamic array?  Odin book says use slices.
	}

	for x in 0..<canv.width{
		for y in 0..<canv.height {
			ray := ray_for_pixel(camera, x, y)
			color := color_at(world, ray)
			set_color(&canv, y, x, color)
		}
	}
	return canv
}

// function render(camera, world)
// image ← canvas(camera.hsize, camera.vsize)
// for y ← 0 to camera.vsize - 1
// for x ← 0 to camera.hsize - 1
// ray ← ray_for_pixel(camera, x, y)
// color ← color_at(world, ray)
// write_pixel(image, x, y, color)
// end for
// end for
// return image
// end function

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
	// TODO: put the over point in comps?  comps.over_point ← comps.point + comps.normalv * EPSILON

	shadowed := is_shadowed(world, comps.point + comps.normalv * EPSILON)
	return lighting(comps.object.material, world.light, comps.point, comps.eyev, comps.normalv, shadowed)
}

color_at :: proc(w: World, r: Ray) -> Color {
	intersections := intersect_world(w, r)
	hit_intersection, ok := hit(intersections)
	if !ok {
		return Color{0,0,0}
	}
	comps := prepare_computations(hit_intersection, r)
	return shade_hit(w, comps)
}

@(test)
test_shade_hit :: proc(t: ^testing.T) {
	r := Ray{make_pnt3(0,0,-5), make_vec3(0,0,1)}
	i := Intersection{4, DefaultWorld.spheres[0]}
	c := shade_hit(DefaultWorld, prepare_computations(i, r))
	testing.expect(t, linalg.vector_length(c - Color{0.38066, 0.47583, 0.2855}) < f64(EPSILON))


	w := DefaultWorld // TODO: Is this a copy? Try in debugger to find out...
	w.light = LightPoint{make_pnt3(0, 0.25, 0), Color{1,1,1}}
	r = Ray{make_pnt3(0, 0, 0), make_vec3(0, 0, 1)}
	i = Intersection{0.5, w.spheres[1]}

	c = shade_hit(w, prepare_computations(i, r))
	testing.expect(t, linalg.vector_length(c - Color{0.90498, 0.90498, 0.90498}) < f64(EPSILON))
}

is_shadowed :: proc(world: World, point: [4]f64) -> bool {
	v := world.light.position - point
	distance := linalg.length(v.xyz)
	direction := linalg.normalize(v)
	r := Ray{point, direction}
	intersections := intersect_world(world, r)
	h, ok := hit(intersections)

	if ok && (h.t < distance) {
		return true
	} else {
		return false
	}
}

lighting :: proc(material: Material, light: LightPoint, point: [4]f64, eyev: [4]f64, normalv: [4]f64, in_shadow: bool) -> Color {
	effective_color := material.color * light.intensity
	lightv := linalg.normalize(light.position - point)
	ambient := effective_color * material.ambient
	if in_shadow {
		return ambient
	}

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

	closest_t := f64(max(f64)) 
	found := false
	hit_record := Intersection{}
	for i in intersections {
		// We only care about intersections that happen in front of the ray (t > 0)
		// and are closer than anything we've found so far.
		if i.t > 0 && i.t < closest_t {
			closest_t = i.t
			hit_record = i
			found = true
		}
	}	
	return hit_record, found
}

position :: proc(ray: Ray, t: f64) -> [4]f64 {
	return ray.origin + ray.direction * t
}

set_color :: proc(canvas: ^Canvas, row: i32, col: i32, color: Color) {
	clipped_color := Color{min(color.r, 1), min(color.g, 1), min(color.b, 1)}
	if (row < canvas.height) && (col < canvas.width) && (row >= 0) && (col >= 0) {
		canvas.pixels[canvas.width*row + col] = clipped_color
	}

}

make_canvas :: proc(width: i32, height: i32) -> Canvas {
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

view_transform :: proc(from: [4]f64, to: [4]f64, up: [4]f64) -> matrix[4,4]f64 {
	// from, to : Pnt3, up : Vec3
	forward := linalg.normalize((to.xyz - from.xyz))
	upn := linalg.normalize(up.xyz)
	left := linalg.cross(forward, upn)
	true_up := linalg.cross(left, forward)
	orientation := matrix[4,4]f64 {
		left.x, left.y, left.z, 0,
		true_up.x, true_up.y, true_up.z, 0,
		-forward.x, -forward.y, -forward.z, 0,
		0, 0, 0, 1
	}
	return orientation * linalg.matrix4_translate([3]f64{-from.x, -from.y, -from.z})
}

main :: proc() {

	floor := make_sphere(0)
	floor.transform = linalg.matrix4_scale([3]f64{10, 0.01, 10})
	floor.material = DefaultMaterial
	floor.material.color = Color{1, 0.9, 0.9}
	floor.material.specular = 0

	left_wall := make_sphere(1)
	left_wall.transform = linalg.matrix4_translate([3]f64{0,0,5}) * linalg.matrix4_rotate(-math.PI/4, [3]f64{0,1,0}) * linalg.matrix4_rotate(math.PI/2, [3]f64{1,0,0}) * linalg.matrix4_scale([3]f64{10, 0.01, 10})

	right_wall := make_sphere(2)
	right_wall.transform = linalg.matrix4_translate([3]f64{0,0,5}) * linalg.matrix4_rotate(math.PI/4, [3]f64{0,1,0}) * linalg.matrix4_rotate(math.PI/2, [3]f64{1,0,0}) * linalg.matrix4_scale([3]f64{10, 0.01, 10})
	right_wall.material = floor.material

	middle := make_sphere(3)
	middle.transform = linalg.matrix4_translate([3]f64{-0.5,1,0.5})
	middle.material = DefaultMaterial
	middle.material.color = Color{0.1, 1, 0.5}
	middle.material.diffuse = 0.7
	middle.material.specular = 0.3

	right := make_sphere(4)
	right.transform = linalg.matrix4_translate([3]f64{1.5, 0.5, -0.5}) * linalg.matrix4_scale([3]f64{0.5,0.5,0.5})
	right.material = DefaultMaterial
	right.material.color = Color{0.5,1,0.1}
	right.material.diffuse = 0.7
	right.material.specular = 0.3

	left := make_sphere(5)
	left.transform = linalg.matrix4_translate([3]f64{-1.5, 0.33, -0.75}) * linalg.matrix4_scale([3]f64{0.33,0.33,0.33})
	left.material = DefaultMaterial
	left.material.color = Color{1,0.8,0.1}
	left.material.diffuse = 0.7
	left.material.specular = 0.3


	world := World{}
	spheres := [6]Sphere{floor, left_wall, right_wall, middle, right, left}
	world.spheres = spheres[:]
	world.light = LightPoint{make_pnt3(-10, 10, -10), Color{1, 1, 1}}
	camera := make_camera(1000, 500, math.PI/3)
	camera.transform = view_transform(make_pnt3(0, 1.5, -5), make_pnt3(0,1,0), make_vec3(0,1,0))
	canv := render(camera, world)

	// canvas_to_ppm(canv, "test.ppm")
	ok := image.write_png("output.png", w=i32(camera.hsize), h=i32(camera.vsize), comp=4, data=raw_data(canvas_to_bmp(canv)), stride_in_bytes=i32(camera.hsize) * 4)

}

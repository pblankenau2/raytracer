package raytracer

import "core:sort"
import "core:fmt"
import "core:os"
import "core:math/linalg"
import "core:math"
import "core:strings"
import "core:slice"
import "core:testing"
import "core:mem/virtual"
import "vendor:stb/image"

EPSILON :: 0.0001

ColorU8 :: [4]u8

canvas_to_bmp :: proc(canv: Canvas) -> [dynamic]ColorU8 {
	pixels := make([dynamic]ColorU8, 0, len(canv.pixels))
	for p, i in canv.pixels {
		append(&pixels, ColorU8{u8(p.r*255), u8(p.g*255), u8(p.b*255), 255})
	}
	return pixels
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
	pixels: []Color
}

Ray :: struct {
	origin: [4]f64, // Point
	direction: [4]f64 // Vec
}

// TODO: Consider making Object struct with a variant union to switch on and contain specialized data
// TODO: Consider just making Object a struct and then defining type aliases (Plane :: Object).  If you need specialized data
// you can later leverage subtyping with the using keyword.  But will switch statements break?
Object :: struct {
	obj_id: i32,
	transform: matrix[4,4]f64,
	material: Material,
	variant: ObjectVariant
}

ObjectVariant :: union {
	Sphere,
	Plane
}
// TODO: consider a SOA instead.
Sphere :: struct {}

// TODO: consider pre-computing the inverse transform and transpose - need to make sure spheres are only created with this function then.
Plane :: struct {}

make_sphere :: proc(
	obj_id: i32,
	transform: matrix[4,4]f64 = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1},
	material : Material = DefaultMaterial
) -> Object {
	// TODO: inv := linalg.inverse(transform)
    // TODO: inv_tr := linalg.transpose(inv)
	return Object{obj_id=obj_id, transform=transform, material=material, variant=Sphere{}}
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
		pixels=make([]Color, camera.hsize*camera.vsize)
	}

	arena: virtual.Arena
	if err := virtual.arena_init_growing(&arena); err != nil {
		fmt.println("Failed to init arena:", err)
		return canv
	}
	defer virtual.arena_destroy(&arena)
	arena_allocator := virtual.arena_allocator(&arena)
	{
		context.allocator = arena_allocator // Set the context allocator for this scope.
		for x in 0..<canv.width{
			for y in 0..<canv.height {
				virtual.arena_free_all(&arena) // Free memory for next pixel
				
				ray := ray_for_pixel(camera, x, y)
				color := color_at(world, ray)
				set_color(canv, y, x, color)
			}
		}
	}
	return canv
}


Intersection :: struct {
	t: f64,
	object: Object, // TODO: make it a generic Object type.
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
	shininess: f64,
	reflective: f64,
	transparency: f64,
	refractive_index: f64,
	pattern: Pattern
}

DefaultMaterial := Material{Color{1,1,1}, 0.1, 0.9, 0.9, 200.0, 0.0, 0.0, 1.0, {}}

World :: struct {
	// TODO. We might want to have struct of arrays here.  Array of spheres, arrays other objects.  Maybe objects should be a tagged union?
	// Also need a count or length field so we can iterate over the arrays that have the same length.
	light: LightPoint,
	objects: []Object
}


DefaultWorld := World{
	LightPoint{make_pnt3(-10,10,-10),Color{1,1,1}},
	{
		make_sphere(
			0, 
			material=Material{
				color=Color{0.8,1.0,0.6},
				ambient=0.1,
				diffuse=0.7,
				specular=0.2,
				shininess=200,
				reflective=0.0,
				transparency=0.0,
				refractive_index=1.0,
				pattern={}
			}
		),
		make_sphere(1, transform=linalg.matrix4_scale([3]f64{0.5,0.5,0.5}))
	}
}


intersect_world :: proc(world: World, ray: Ray) -> [dynamic]Intersection {
	intersections := make([dynamic]Intersection, 0, 32) // TODO: Preallocate 2*number of objects.
	for s, i in world.objects {
		i1, i2, okay := intersect(ray, s)
		if okay {
			append(&intersections, i1, i2)
		}
	}
	slice.sort_by(intersections[:], proc(i, j: Intersection) -> bool {
        return i.t < j.t
    })

    return intersections
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
	object: Object,
	point: [4]f64,
	eyev: [4]f64,
	normalv: [4]f64,
	inside: bool,
	reflectv: [4]f64,
	n1: f64, // refractive_index_exiting
	n2: f64, // refractive_index_entering
	under_point: [4]f64,
}

prepare_computations :: proc(intersection: Intersection, ray: Ray, xs: []Intersection) -> PreComputations {
	point := position(ray, intersection.t)
	normalv := normal_at(intersection.object, point)
	eyev := -ray.direction
	reflectv := reflect(ray.direction, normalv)

	inside := false
	if linalg.dot(normalv, eyev) < 0 {
		inside = true
		normalv = -normalv
	}

	containers : [dynamic]Intersection
	hit_x, ok := hit(xs)
	n1 := 1.0
	n2 := 1.0
	for x in xs {
		if x == hit_x {
			if (len(containers) == 0) {
				n1 = 1.0
			} else {
				n1 = containers[len(containers) - 1].object.material.refractive_index
			}
		}

		index, found := slice.linear_search(containers[:], x)
		if found {
			// containers[index] = Intersection{}
			ordered_remove(&containers, index)
		} else {
			append(&containers, x)
		}
		if x == hit_x {
			if len(containers) == 0 {
				n2 = 1.0
			} else {
				n2 = containers[len(containers) - 1].object.material.refractive_index
			}
		}

	}

	return PreComputations{
		t=intersection.t,
		object=intersection.object,
		point=point,
		eyev=eyev,
		normalv=normalv,
		inside=inside,
		reflectv=reflectv,
		n1=n1,
		n2=n2,
		under_point=point - normalv * EPSILON
	}
}

shade_hit :: proc(world: World, comps: PreComputations, remaining: int) -> Color {
	// TODO: put the over point in comps?  comps.over_point ← comps.point + comps.normalv * EPSILON
	over_point := comps.point + comps.normalv * EPSILON
	shadowed := is_shadowed(world, over_point)
	surface := lighting(comps.object.material, comps.object.transform, world.light, comps.point, comps.eyev, comps.normalv, shadowed)

	reflect_ray := Ray{over_point, comps.reflectv}
	if remaining <= 0 {
		return surface
	}
	color := color_at(world, reflect_ray, remaining-1)
	reflected := color * comps.object.material.reflective

	// TODO: call refracted color here.

	return surface + reflected
}

color_at :: proc(w: World, r: Ray, remaining: int = 5) -> Color {
	intersections := intersect_world(w, r)
	hit_intersection, ok := hit(intersections[:])
	if !ok {
		return Color{0,0,0}
	}
	comps := prepare_computations(hit_intersection, r, intersections[:])
	return shade_hit(w, comps, remaining)
}

@(test)
test_shade_hit :: proc(t: ^testing.T) {
	r := Ray{make_pnt3(0,0,-5), make_vec3(0,0,1)}
	i := Intersection{4, DefaultWorld.objects[0]}
	c := shade_hit(DefaultWorld, prepare_computations(i, r, []Intersection{i}), 5)
	testing.expect(t, linalg.vector_length(c - Color{0.38066, 0.47583, 0.2855}) < f64(EPSILON))


	w := DefaultWorld
	w.light = LightPoint{make_pnt3(0, 0.25, 0), Color{1,1,1}}
	r = Ray{make_pnt3(0, 0, 0), make_vec3(0, 0, 1)}
	i = Intersection{0.5, w.objects[1]}

	c = shade_hit(w, prepare_computations(i, r, []Intersection{i}), 5)
	testing.expect(t, linalg.vector_length(c - Color{0.90498, 0.90498, 0.90498}) < f64(EPSILON))
}

is_shadowed :: proc(world: World, point: [4]f64) -> bool {
	v := world.light.position - point
	distance := linalg.length(v.xyz)
	direction := linalg.normalize(v)
	r := Ray{point, direction}
	intersections_shadowed := intersect_world(world, r)
	h, ok := hit(intersections_shadowed[:])

	if ok && (h.t < distance) {
		return true
	} else {
		return false
	}
}

Stripe :: struct {
	color1: Color,
	color2: Color,
	transform: matrix[4,4]f64
}

Gradient :: struct {
	color1: Color,
	color2: Color,
	transform: matrix[4,4]f64
}

Checker :: struct {
	color1: Color,
	color2: Color,
	transform: matrix[4,4]f64
}

Pattern :: union {
	Stripe,
	Gradient,
	Checker,
}

lighting :: proc(material: Material, object_transform: matrix[4,4]f64, light: LightPoint, point: [4]f64, eyev: [4]f64, normalv: [4]f64, in_shadow: bool) -> Color {
	effective_color := material.color * light.intensity
	object_point := linalg.inverse(object_transform) * point
	switch pattern in material.pattern {
		case Stripe:
			pattern_point := linalg.inverse(pattern.transform) * object_point
			if (int(math.floor(pattern_point.x)) %% 2) == 1 {
				effective_color = pattern.color1 * light.intensity
			} else {
				effective_color = pattern.color2 * light.intensity
			}
		case Checker:
			pattern_point := linalg.inverse(pattern.transform) * object_point
			if (int(math.floor(pattern_point.x) + math.floor(pattern_point.y) + math.floor(pattern_point.z)) %% 2) == 1 {
				effective_color = pattern.color1 * light.intensity
			} else {
				effective_color = pattern.color2 * light.intensity
			}
		case Gradient:
			pattern_point := linalg.inverse(pattern.transform) * object_point

			distance := pattern.color2 - pattern.color1
			fraction := pattern_point.x - math.floor(pattern_point.x)
			effective_color = (pattern.color1 + distance * fraction) * light.intensity
	}


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
intersect :: proc(ray: Ray, object: Object) -> (Intersection, Intersection, bool) {
	switch o in object.variant {
	case Sphere:
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
	case Plane:
		if math.abs(ray.direction.y) < EPSILON {
			return {}, {}, false
		}
		t := -ray.origin.y / ray.direction.y
		return Intersection{t, object}, Intersection{t, object}, true
	}
	return {}, {}, false
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

set_color :: proc(canvas: Canvas, row: i32, col: i32, color: Color) {
	clipped_color := Color{min(color.r, 1), min(color.g, 1), min(color.b, 1)}
	if (row < canvas.height) && (col < canvas.width) && (row >= 0) && (col >= 0) {
		// canvas is an immutable reference, so you can't overwrite the pixels slice, but you can right to the 
		// the heap allocated data the slice points to.
		canvas.pixels[canvas.width*row + col] = clipped_color 
	}

}

make_canvas :: proc(width: i32, height: i32) -> Canvas {
	pix := make([]Color, width*height)

	for i in 0..<(width*height) {
		pix[i] = Color{0,0,0}
	}

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

normal_at :: proc(obj: Object, p: [4]f64) -> [4]f64 {

	switch o in obj.variant {
	case Sphere:
		object_point := linalg.inverse(obj.transform) * p
		object_normal := object_point - [4]f64{0.0,0.0,0.0,0.0}
		world_normal := linalg.transpose(linalg.inverse(obj.transform)) * object_normal
		world_normal.w = 0.0
		return linalg.normalize(world_normal)
	case Plane:
		return make_vec3(0, 1, 0)
	}
	return {} // Weird that we need to do this since every case returns.
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

	// floor := make_sphere(0)
	// floor.transform = linalg.matrix4_scale([3]f64{10, 0.01, 10})
	// floor.material = DefaultMaterial
	// floor.material.color = Color{1, 0.9, 0.9}
	// floor.material.specular = 0

	floor := Object{0, linalg.identity_matrix(matrix[4,4]f64), DefaultMaterial, Plane{}}
	floor.material.color = Color{1, 0.9, 0.9}
	floor.material.specular = 0
	floor.material.pattern = Checker{
		Color{0.1, 0.1, 0.1},
		Color{0.99, 0.99, 0.99},
		linalg.identity_matrix(matrix[4,4]f64) * linalg.matrix4_translate([3]f64{0,-EPSILON,0}) // * linalg.matrix4_translate([3]f64{-1.0, 0, 0}) * linalg.matrix4_scale([3]f64{2.0,2.0,2.0})
	}
	floor.material.reflective = 0.5

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
	middle.material.pattern = Gradient{
		Color{0.9, 0.1, 0.1},
		Color{0.1, 1, 0.5},
		// linalg.matrix4_scale([3]f64{0.25,0.25,0.25}) * matrix[4,4]f64{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1}
		linalg.identity_matrix(matrix[4,4]f64) * linalg.matrix4_translate([3]f64{-1.0, 0, 0}) * linalg.matrix4_scale([3]f64{2.0,2.0,2.0})
	}

	right := make_sphere(4)
	right.transform = linalg.matrix4_translate([3]f64{1.5, 0.5, -0.5}) * linalg.matrix4_scale([3]f64{0.5,0.5,0.5})
	right.material = DefaultMaterial
	right.material.color = Color{0.8, 0.8, 0.8}
	right.material.diffuse = 0.3
	right.material.specular = 0.7
	right.material.shininess = 300
	right.material.reflective = 0.95

	left := make_sphere(5)
	left.transform = linalg.matrix4_translate([3]f64{-1.5, 0.33, -0.75}) * linalg.matrix4_scale([3]f64{0.33,0.33,0.33})
	left.material = DefaultMaterial
	left.material.color = Color{0.6352941176470588, 0.5098039215686274, 0.7607843137254902}
	left.material.diffuse = 0.7
	left.material.specular = 0.3


	world := World{}

	// objects := [6]Object{floor, left_wall, right_wall, middle, right, left}
	objects := [4]Object{floor, middle, right, left}
	world.objects = objects[:]
	world.light = LightPoint{make_pnt3(-10, 10, -10), Color{1, 1, 1}}
	camera := make_camera(1000, 500, math.PI/3)
	camera.transform = view_transform(make_pnt3(0, 1.5, -5), make_pnt3(0,1,0), make_vec3(0,1,0))
	canv := render(camera, world)

	// canvas_to_ppm(canv, "test.ppm")
	ok := image.write_png("output.png", w=i32(camera.hsize), h=i32(camera.vsize), comp=4, data=raw_data(canvas_to_bmp(canv)), stride_in_bytes=i32(camera.hsize) * 4)

}

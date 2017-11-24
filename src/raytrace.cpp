#include "scene.h"

#include <thread>
#include "ext/yocto_utils.h"


#ifndef FALSE
	#define FALSE 0
#endif
#ifndef TRUE
	#define TRUE 1
#endif


#if !defined(DEBUGBREAK)
	#if defined(__clang__)
		#define DEBUGBREAK __builtin_trap()
	#elif defined(__GNUC__)
		#define DEBUGBREAK __asm__ __volatile__ ("int3")
	#elif defined(_MSC_VER)
		#define DEBUGBREAK __debugbreak()
	#endif
#endif

#if ENABLE_ASSERT
	#define ASSERT(c) do{ if(!(c)) DEBUGBREAK; } while(0)
#else
	#define ASSERT(c) do{ (void)sizeof(c); } while(0)
#endif


#if !DISABLE_MULTITHREADING
	#include <thread>
#endif


static
vec3f
kajiya_kay(const scene* scn, vec3f ws_point, vec3f normal, vec3f ws_camera_pos,
           vec3f light_intensity, vec3f ws_light_pos,
           vec3f diffuse_color, vec3f specular_color, float specular_exp)
{
	constexpr float RAY_EPSILON = 1e-2f;
	
	vec3f ps_light_pos   = ws_light_pos - ws_point;
	vec3f ps_light_dir   = normalize(ps_light_pos);
	float light_distance = length(ps_light_pos);
	
#if DISABLE_TRANSPARENCY
	constexpr float transparency = 1.0f;
	ray3f ray = {ws_point, ps_light_dir, RAY_EPSILON, light_distance-RAY_EPSILON};
	bool blocked = intersect_any(scn, ray);
	if(blocked)
		return {0.0f, 0.0f, 0.0f};
#else
	float transparency = 1.0f;
	ray3f ray = {ws_point, ps_light_dir, RAY_EPSILON, light_distance-RAY_EPSILON};
	intersection3f intersection;
	do
	{
		intersection = intersect_first(scn, ray);
		if(intersection.hit())
		{		
			instance* hit_object = intersection.ist;
			if(hit_object->mat)
			{
				transparency *= (1.0f - hit_object->mat->op);
				if(transparency <= 0.0f)
					return {0.0f, 0.0f, 0.0f};
			}
			ray.o    += intersection.dist * ray.d;
			ray.tmax -= intersection.dist;
		}
	} while(intersection.hit());
#endif
	
	vec3f ps_camera_pos = ws_camera_pos - ws_point;
	vec3f ps_camera_dir = normalize(ps_camera_pos);
	vec3f v_l_bisector  = normalize(ps_camera_dir + ps_light_dir);
	
	float coeff_n_l = fabsf(dot(normal, ps_light_dir));
	float coeff_n_h = fabsf(dot(normal, v_l_bisector));
	float diff_angle_coefficient = sqrtf(1.0f - coeff_n_l*coeff_n_l);
	float spec_angle_coefficient = sqrtf(1.0f - coeff_n_h*coeff_n_h);
	
	float kt = transparency;
	vec3f kd = diffuse_color;
	vec3f ks = specular_color;
	vec3f I  = light_intensity;
	float r2 = light_distance*light_distance;
	
	vec3f color = {};
	color += kt * kd * (I / r2) *      max(0.0f, diff_angle_coefficient);
	color += kt * ks * (I / r2) * powf(max(0.0f, spec_angle_coefficient), specular_exp);
	return color;
}


static
vec3f
blinn_phong(const scene* scn, vec3f ws_point, vec3f normal, vec3f ws_camera_pos,
            vec3f light_intensity, vec3f ws_light_pos,
            vec3f diffuse_color, vec3f specular_color, float specular_exp)
{
	constexpr float RAY_EPSILON = 1e-4f;
	
	vec3f ps_light_pos   = ws_light_pos - ws_point;
	vec3f ps_light_dir   = normalize(ps_light_pos);
	float light_distance = length(ps_light_pos);

#if DISABLE_TRANSPARENCY
	constexpr float transparency = 1.0f;
	ray3f ray = {ws_point, ps_light_dir, RAY_EPSILON, light_distance-RAY_EPSILON};
	bool blocked = intersect_any(scn, ray);
	if(blocked)
		return {0.0f, 0.0f, 0.0f};
#else
	float transparency = 1.0f;
	ray3f ray = {ws_point, ps_light_dir, RAY_EPSILON, light_distance-RAY_EPSILON};
	intersection3f intersection;
	do
	{
		intersection = intersect_first(scn, ray);
		if(intersection.hit())
		{		
			instance* hit_object = intersection.ist;
			if(hit_object->mat)
			{
				transparency *= (1.0f - hit_object->mat->op);
				if(transparency <= 0.0f)
					return {0.0f, 0.0f, 0.0f};
			}
			ray.o    += intersection.dist * ray.d;
			ray.tmax -= intersection.dist;
		}
	} while(intersection.hit());
#endif
	
	vec3f ps_camera_pos = ws_camera_pos - ws_point;
	vec3f ps_camera_dir = normalize(ps_camera_pos);
	vec3f v_l_bisector  = normalize(ps_camera_dir + ps_light_dir);
	
	float diff_angle_coefficient = dot(normal, ps_light_dir);
	float spec_angle_coefficient = dot(normal, v_l_bisector);
	
	float kt = transparency;
	vec3f kd = diffuse_color;
	vec3f ks = specular_color;
	vec3f I  = light_intensity;
	float r2 = light_distance*light_distance;
	
	vec3f color = {};
	color += kt * kd * (I / r2) *      max(0.0f, diff_angle_coefficient);
	color += kt * ks * (I / r2) * powf(max(0.0f, spec_angle_coefficient), specular_exp);
	return color;
}


ray3f
eval_camera(const camera* cam, const vec2f& uv)
{
	constexpr float d = 1.0f; // Always fixed d?
	
	float h = d * 2.0f * tanf(cam->fovy*0.5f);
	float w = d * h    * cam->aspect;
	vec3f direction = (uv.u-0.5f) * w * cam->frame.x +
	                  (uv.v-0.5f) * h * cam->frame.y -
	                  d * cam->frame.z;
	ray3f ray = {cam->frame.o, normalize(direction)};
	return ray;
}


vec3f
lookup_texture(const texture* txt, int x, int y, bool srgb)
{
	ASSERT(txt);
	ASSERT(txt->ldr.pixels.size() > 0);
	ASSERT(x < txt->ldr.width);
	ASSERT(y < txt->ldr.height);
	
	constexpr float INV_255 = 1.0f / 255.0f;
	constexpr float GAMMA   = 2.2f;
	
	vec4b texel_255 = txt->ldr.pixels[y*txt->ldr.width + x];
	vec4f texel     = {texel_255.r * INV_255,
		               texel_255.g * INV_255,
		               texel_255.b * INV_255,
		               texel_255.a * INV_255};
	if(srgb)
	{
		texel.r = powf(texel.r, GAMMA);
		texel.g = powf(texel.g, GAMMA);
		texel.b = powf(texel.b, GAMMA);
	}
	
	// TODO: Check this
	// Premultiply
	vec4f result = texel * texel.a;
	return result.rgb;
}


vec3f
eval_texture(const texture* txt, const vec2f& texcoord, bool srgb)
{
	float ignored;
	float v = modff(texcoord.v, &ignored) * txt->ldr.height;
	float u = modff(texcoord.u, &ignored) * txt->ldr.width;
	
	int texel_x0 = floor(u);
	int texel_y0 = floor(v);
#if 0
	int texel_x0 = (int)u;
	int texel_y0 = (int)v;
#endif
	int texel_x1 = (texel_x0 + 1) % txt->ldr.width;
	int texel_y1 = (texel_y0 + 1) % txt->ldr.height;
	
	vec3f sample_x0y0 = lookup_texture(txt, texel_x0, texel_y0, srgb);
	vec3f sample_x1y0 = lookup_texture(txt, texel_x1, texel_y0, srgb);
	vec3f sample_x0y1 = lookup_texture(txt, texel_x0, texel_y1, srgb);
	vec3f sample_x1y1 = lookup_texture(txt, texel_x1, texel_y1, srgb);
	
	float wx = u - texel_x0;
	float wy = v - texel_y0;
	
	sample_x0y0 *= (1.0f-wx) * (1.0f-wy);
	sample_x1y0 *= (     wx) * (1.0f-wy);
	sample_x0y1 *= (1.0f-wx) * (     wy);
	sample_x1y1 *= (     wx) * (     wy);
	
	vec3f color = sample_x0y0 + sample_x1y0 + sample_x0y1 + sample_x1y1;
	return color;
}


vec4f
shade(const scene* scn, const std::vector<instance*>& lights, const vec3f& ambient_light, const ray3f& ray)
{
	vec4f result = {0.0f, 0.0f, 0.0f, 1.0f};
	
	intersection3f intersection = intersect_first(scn, ray);
	if(intersection.hit())
	{
		ASSERT(intersection.ist);
		ASSERT(intersection.ist->shp);
		
		instance* hit_object          = intersection.ist;
		int       hit_primitive_index = intersection.ei;
		vec4f     hit_barycentric     = intersection.ew;
		
		// ASSERT(hit_object->name == "floor");
		if(hit_object->mat)
		{
			material* material = hit_object->mat;
			
			vec3f hit_position = eval_pos(hit_object, hit_primitive_index, hit_barycentric);
			vec3f hit_normal   = eval_norm(hit_object, hit_primitive_index, hit_barycentric);
			vec2f hit_texcoord = eval_texcoord(hit_object->shp, hit_primitive_index, hit_barycentric);
			
			vec3f diff_color = material->kd_txt ? eval_texture(material->kd_txt, hit_texcoord, TRUE) : material->kd;
			vec3f spec_color = material->ks_txt ? eval_texture(material->ks_txt, hit_texcoord, TRUE) : material->ks;
			vec3f refl_color = material->kr_txt ? eval_texture(material->kr_txt, hit_texcoord, TRUE) : material->kr;
			
			// Apply ambient
			vec3f color = diff_color * ambient_light;
			
			// Apply diffuse + specular
			vec3f camera_pos   = ray.o;
			float specular_exp = material->rs ? 2.0f / powf(material->rs, 4.0f) - 2.0f : 1e6f;
			
			if(hit_object->shp->lines.size() > 0)
			{
				for(instance* light: lights)
				{
					ASSERT(light->mat && light->shp && light->shp->pos.size());
					
					vec3f light_intensity = light->mat->ke;
					vec3f light_position  = transform_point(light->frame, light->shp->pos[0]);
					color += kajiya_kay(scn, hit_position, hit_normal, camera_pos,
					                    light_intensity, light_position,
					                    diff_color, spec_color, specular_exp);
				}
			}
			else // triangles
			{
				for(instance* light: lights)
				{
					ASSERT(light->mat && light->shp && light->shp->pos.size());
					
					vec3f light_intensity = light->mat->ke;
					vec3f light_position  = transform_point(light->frame, light->shp->pos[0]);
					color += blinn_phong(scn, hit_position, hit_normal, camera_pos,
					                     light_intensity, light_position,
					                     diff_color, spec_color, specular_exp);
				}
			}
			
			// Apply reflection
			if(refl_color.r > 0.0f
			|| refl_color.g > 0.0f
			|| refl_color.b > 0.0f)
			{	
				vec3f ps_camera_dir   = normalize(camera_pos - hit_position);
				vec3f reflected_dir   = 2.0f * dot(hit_normal, ps_camera_dir) * hit_normal - ps_camera_dir;
				ray3f reflected_ray   = {hit_position, reflected_dir};
				vec4f reflected_color = shade(scn, lights, ambient_light, reflected_ray);
				color += refl_color * reflected_color.rgb;
			}
			
#if !DISABLE_TRANSPARENCY
			float opacity = material->op;
			if(opacity < 1.0f)
			{
				ray3f ray_beyond   = {hit_position, ray.d};
				vec4f color_beyond = shade(scn, lights, ambient_light, ray_beyond);
				color = opacity * color + (1.0f - opacity) * color_beyond.rgb;
			}
#endif
			
			result.rgb = color;
		}
		else // !hit_object->mat
		{			
			result.rgb = {1.0f, 0.0f, 1.0f}; // Bright purple - everyone's favorite error color
		}
	}
	
	return result;
}


void
raytrace_region(const scene* scn, vec3f ambient_light, int samples,
                const camera* cam, image4f* image, const std::vector<instance*>& lights,
                int x0, int y0, int x1, int y1)
{
	float w = image->width;
	float h = image->height;
	for(int y = y0; y < y1; ++y)
	{
		for(int x = x0; x < x1; ++x)
		{
			vec4f color = {};
			for(int sy = 0; sy < samples; ++sy)
			{
				for(int sx = 0; sx < samples; ++sx)
				{
					vec2f uv     = {       (x + (sx+0.5f) / samples) / w,
						            1.0f - (y + (sy+0.5f) / samples) / h};
					ray3f ray    = eval_camera(cam, uv);
					vec4f sample = shade(scn, lights, ambient_light, ray);
					
					color += sample;
				}
			}
			image->pixels[y*image->width + x] = color / (samples*samples);
		}
	}
}


image4f
raytrace(const scene* scn, const vec3f& amb, int resolution, int samples)
{
	auto cam = scn->cameras.front();
	auto img = image4f((int)std::round(cam->aspect * resolution), resolution);
	
#if DISABLE_MULTISAMPLING
	samples = 1;
#endif

	// IL TUO CODICE VA QUI
	std::vector<instance*> lights = {};
	for(instance* object: scn->instances)
	{
		if(object->mat)
		{
			vec3f emission = object->mat->ke;
			if(emission.r > 0.0f
			|| emission.g > 0.0f
			|| emission.b > 0.0f)
			{
				lights.push_back(object);
			}
		}
	}
	
#if DISABLE_MULTITHREADING
	raytrace_region(scn, amb, samples,
	                cam, &img, lights,
	                0, 0, img.width, img.height);
#else
	std::thread* threads       = NULL;
	unsigned int threads_count = std::thread::hardware_concurrency();
	if(threads_count > 1)
		threads = new std::thread[threads_count-1];
	
	int height_per_thread = img.height / threads_count;
	int leftover_height   = img.height % threads_count;
	
	for(unsigned int ti = 1; ti < threads_count; ++ti)
	{
		int y0 = ti * height_per_thread;
		int y1 = y0 + height_per_thread;
		if(leftover_height > 0)
			++y1, --leftover_height;
		
		threads[ti-1] = std::thread(raytrace_region,
		                            scn, amb, samples,
		                            cam, &img, std::ref(lights), 
		                            0, y0, img.width, y1);
	}
	
	ASSERT(leftover_height == 0);
	raytrace_region(scn, amb, samples,
	                cam, &img, lights,
	                0, 0, img.width, height_per_thread);
	

	for(unsigned int ti = 1; ti < threads_count; ++ti)
	{
		threads[ti-1].join();
	}
	
	if(threads)
		delete[] threads;
#endif

	return img;
}


int
main(int argc, char** argv) {
	// command line parsing
	auto parser =
		yu::cmdline::make_parser(argc, argv, "raytrace", "raytrace scene");
	auto resolution = yu::cmdline::parse_opti(
		parser, "--resolution", "-r", "vertical resolution", 720);
	auto samples = yu::cmdline::parse_opti(
		parser, "--samples", "-s", "per-pixel samples", 1);
	auto amb = yu::cmdline::parse_optf(
		parser, "--ambient", "-a", "ambient color", 0.1f);
	auto imageout = yu::cmdline::parse_opts(
		parser, "--output", "-o", "output image", "out.png");
	auto scenein = yu::cmdline::parse_args(
		parser, "scenein", "input scene", "scene.obj", true);
	yu::cmdline::check_parser(parser);

	// load scene
	printf("loading scene %s\n", scenein.c_str());
	auto scn = load_scene(scenein);

	// create bvh
	printf("creating bvh\n");
	build_bvh(scn, false);

	// raytrace
	printf("tracing scene\n");
	auto hdr = raytrace(scn, vec3f{amb, amb, amb}, resolution, samples);

	// tonemap and save
	printf("saving image %s\n", imageout.c_str());
	save_hdr_or_ldr(imageout, hdr);
}

  #pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!

    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);

    if (!vertex_) {
        // Hit background. Account for the environment map if needed.
        if (has_envmap(scene)) {
            const Light& envmap = get_envmap(scene);
            return emission(envmap,
                -ray.dir, // pointing outwards from light
                ray_diff.spread,
                PointAndNormal{}, // dummy parameter for envmap
                scene);
        }
        return make_zero_spectrum();
    }

    PathVertex vertex = *vertex_;
    Spectrum radiance = make_zero_spectrum();
	Spectrum sigma_a = get_sigma_a(scene.media[vertex.exterior_medium_id], vertex.position);
	double t = distance(vertex.position, ray.org);
	Spectrum transmittance = exp(-sigma_a * t);

	Spectrum Le(0, 0, 0);
    if(is_light(scene.shapes[vertex.shape_id]))
        Le = emission(vertex, -ray.dir, scene);

    return transmittance * Le;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);


	Real u = next_pcg32_real<Real>(rng);
	Spectrum sigma_a = get_sigma_a(scene.media[scene.camera.medium_id], ray.org);
	Spectrum sigma_s = get_sigma_s(scene.media[scene.camera.medium_id], ray.org);
	Spectrum sigma_t = sigma_a + sigma_s;
	// sample t s.t. p(t) ~ exp(-sigma_t * t)
	Real t = log(1 - u) / (-sigma_t[0]);
    // if ray doesn't hit surface geometry
    Real t_hit = infinity<Real>();
    if(vertex_)
        t_hit = distance(vertex_->position, ray.org);

    // Res
	Spectrum L_2 = make_zero_spectrum();

	// eqn (9) has 2 parts, sample t to decide which part we are in.
    if (t < t_hit) {
        // pdf eqn (11)
		Spectrum pdf_transmittance = sigma_t * exp(-sigma_t * t);
		Spectrum transmittance = exp(-sigma_t * t);     // transmittance from ray origin to scatter point
        
		// compute the point where scattering happens
		Vector3 scatter_point = ray.org + ray.dir * t;
        // sampling light source
        Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
        Real light_w = next_pcg32_real<Real>(rng);
        Real shape_w = next_pcg32_real<Real>(rng);
        int light_id = sample_light(scene, light_w);
        const Light& light = scene.lights[light_id];
        PointAndNormal point_on_light =
            sample_point_on_light(light, scatter_point, light_uv, shape_w, scene);

		Vector3 wi = normalize(point_on_light.position - scatter_point);
		Real dis_sp_l = distance(scatter_point, point_on_light.position);   // scatter point to light
		Real G = (abs(dot(-wi, point_on_light.normal))) / (dis_sp_l * dis_sp_l);
		PhaseFunction pf = get_phase_function(scene.media[scene.camera.medium_id]);
		Spectrum rho = eval(pf, wi, -ray.dir);
		Spectrum trans_pp = exp(-sigma_t * dis_sp_l);   // transmittance from scatter point to light
		Vector3 dir_sp_light = normalize(point_on_light.position - scatter_point);
        Ray shadow_ray{ scatter_point, dir_sp_light,
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, scatter_point) };
        Real V = 0;
        if (!occluded(scene, shadow_ray))
            V = 1;
        
		Spectrum Le = emission(light, -dir_sp_light, Real(0), point_on_light, scene);
        Real pdf_light = light_pmf(scene, light_id) *
			pdf_point_on_light(light, point_on_light, scatter_point, scene);
        // monte carlo for L_scatter
		Spectrum L_scatter1 = rho * Le * trans_pp * G * V / pdf_light;      

		// monte carlo for L_2 in eqn (9)
		L_2 += transmittance * sigma_s * L_scatter1 / pdf_transmittance;
    }
    // sampled t > t_hit, we treat it as hitting the surface.
    else {
		// monte carlo for hitting surface contribution
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum pdf_trans = exp(-sigma_t * t_hit);
		assert(vertex_.has_value());
        if (is_light(scene.shapes[vertex_->shape_id])) {
			L_2 += (transmittance * emission(*vertex_, -ray.dir, scene)) / pdf_trans;
        }
        
    }
    return L_2;
}

// update medium id the ray currently lie in
int update_medium_id(const PathVertex& isect, const Ray& ray, int medium_id) {
    if (isect.exterior_medium_id != isect.interior_medium_id) {
        // At medium transition. Update Medium
        if (dot(ray.dir, isect.geometric_normal) > 0) {
            return isect.exterior_medium_id;
        }
        else
            return isect.interior_medium_id;
    }
    return medium_id;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

	int current_medium_id = scene.camera.medium_id;
	Spectrum tp = Vector3{ 1, 1, 1 };   // path throughput
	Spectrum radiance = make_zero_spectrum();
    int depth = 0;

    while (true) {
        bool scatter = false;
		Spectrum transmittance(1, 1, 1);
        Spectrum pdf_trans(1, 1, 1);

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
		Real t_hit = infinity<Real>();
        if (vertex_)
            t_hit = distance(vertex_->position, ray.org);
		Real t = infinity<Real>();


        if (current_medium_id != -1) {
            const Medium& current_medium = scene.media[current_medium_id];
			Spectrum sigma_a = get_sigma_a(current_medium, ray.org);
			Spectrum sigma_s = get_sigma_s(current_medium, ray.org);
			Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);

            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            t = log(1 - u) / (-sigma_t[0]);

			// scatter in the medium
            if (t < t_hit) {
                scatter = true;
                // pdf eqn (11)
                pdf_trans = sigma_t * exp(-sigma_t * t);
                transmittance = exp(-sigma_t * t);     // transmittance from ray origin to scatter point
				ray.org = ray.org + ray.dir * t;
            }
			// hit surface 
            else {
                scatter = false;
                pdf_trans = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
				// avoid self intersection
				assert(vertex_.has_value());
                ray = Ray{ vertex_->position, ray.dir, get_intersection_epsilon(scene), infinity<Real>() };
            }
        }
		tp = tp * (transmittance/pdf_trans);

		// latter part of eqn(14), contribution from hitting emissive geometry surface
        if (!scatter) {
			Spectrum Le(0, 0, 0);

            assert(vertex_.has_value());
            if(is_light(scene.shapes[vertex_->shape_id]))
				Le = emission(*vertex_, -ray.dir, scene);

            radiance += tp * Le;
        }

        // reaching max depth
        if (depth == scene.options.max_depth - 1 && scene.options.max_depth != -1)
            break;

        // if hit surface
        if (!scatter && vertex_) {
            // if the surface is not geometric: hitting another participating medium.
            if (vertex_->material_id == -1) {
                // index matching interface, skip through it
				current_medium_id = update_medium_id(*vertex_, ray, current_medium_id);
                // added when implementing version4
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                depth++;
                continue;
            }
        }

        // sample next dir & update path tp
        if (scatter) {
			PhaseFunction pf = get_phase_function(scene.media[current_medium_id]);
			Vector2 pf_rnd_param{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            std::optional<Vector3> nxt_dir = sample_phase_function(pf, -ray.dir, pf_rnd_param);
			assert(nxt_dir.has_value());
			Real dir_pdf = pdf_sample_phase(pf, -ray.dir, *nxt_dir);

			Spectrum rho = eval(pf, *nxt_dir, -ray.dir);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
			tp = tp * (rho * sigma_s / dir_pdf);
			ray.dir = *nxt_dir;
        }
        else {  // hit surface, do nothing for now
            break;
        }

		// russian roulette
        Real rr_prob = 1;
        if (depth >= scene.options.rr_depth) {
			rr_prob = min(tp[0], 0.95);
            if(next_pcg32_real<Real>(rng) > rr_prob) 
                break;
            else 
				tp = tp / rr_prob;
        }
        depth++;
    }
    return radiance;
}

// next event estimation for participating medium, it is used in version 4 and later.
// wi: output ray, pointing towards the camera path 
Spectrum nee_medium(const Scene& scene, const Vector3 &point, int current_medium_id, pcg32_state& rng, RayDifferential& ray_diff, int bounce, const Vector3& wi) {
	Vector3 p = point;  // copy for later update while walking the ray
    // sample light
    Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light& light = scene.lights[light_id];
    PointAndNormal p_prime =
        sample_point_on_light(light, p, light_uv, shape_w, scene);

    // Compute transmittance to light. Skip through index-matching shapes.
    Spectrum T_light(1, 1, 1);
    int shadow_bounce = 0;
    Real pdf_trans_dir = 1;     // for MIS
	int shadow_medium_id = current_medium_id;
    Real pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, p, scene);

    // this loop is only to update the transmittance && pdf of the path from shading point to the sampled light point.
    while (1) {
		Ray shadow_ray = Ray{ p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime.position, p) };

        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
		Real next_t = distance(p, p_prime.position);
        if (isect) {
			next_t = distance(p, isect->position);
        }

		// account for transmittance to next_t
        if (shadow_medium_id != -1) {
            const Medium& shadow_medium = scene.media[shadow_medium_id];
			Spectrum sigma_a = get_sigma_a(shadow_medium, shadow_ray.org);
			Spectrum sigma_s = get_sigma_s(shadow_medium, shadow_ray.org);
			Spectrum sigma_t = sigma_a + sigma_s;
			T_light *= exp(-sigma_t * next_t);
			pdf_trans_dir *= exp(-sigma_t[0] * next_t);
        }

        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
			// something is blocking, is it opaque surface ?
            if(isect->material_id != -1) {
                // hit opaque surface, blocked
                return make_zero_spectrum();
			}

			// hit index-matching surface, we want to pass through -- this introduces one extra connection vertex
            shadow_bounce++;
            if(scene.options.max_depth != -1 && bounce +  shadow_bounce + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
			}

            shadow_medium_id = update_medium_id(*isect, shadow_ray, shadow_medium_id);
			p = p + next_t * shadow_ray.dir;
        }
    }
    if (T_light[0] > 0) {
		Real G = (abs(dot(-normalize(p_prime.position - point), p_prime.normal))) / (distance_squared(point, p_prime.position));
		Vector3 dir_light = normalize(p_prime.position - point);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);

		PhaseFunction pf = get_phase_function(scene.media[current_medium_id]);
		Spectrum rho = eval(pf, wi, dir_light);
        Spectrum contrib = T_light * rho * G * L / pdf_light;
        // MIS: it’s also possible
        // that a phase function sampling + multiple exponential sampling  will reach the light source.
        // We also need to multiply with G to convert phase function PDF to area measure.
        Real pdf_phase = pdf_sample_phase(pf, wi, dir_light) * G * pdf_trans_dir;

		Real w = (pdf_light * pdf_light) / (pdf_light * pdf_light + pdf_phase * pdf_phase);
		return w * contrib;
    }
	return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    // MIS is only for direct illumination
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    int current_medium_id = scene.camera.medium_id;
    Spectrum tp = Vector3{ 1, 1, 1 };   // path throughput
    Spectrum radiance = make_zero_spectrum();
    int depth = 0;

    /*
        For pdf_nee and pdf_phase, when caculating the MIS weight, the only 2 things I care about are,
        the probability of sampling this light point (getting here) from the last scatter point using phase function sampling:
            I need probability of sampling this dir from phase function, and prob. of sampling the ray length t.
        And
        the probability of sampling this light_point (getting here) from the last scatter point using nee:
			I need probability of sampling this light point, which is pdf_light
    */
    Real pdf_dir_latest = 0;           // In solid angle measure: the pdf of the latest phase function sampling
    Real pdf_multi_trans = 1;          // The product PDF of transmittance sampling going through several index-matching surfaces
                                       // from the last phase function sampling
	Vector3 pos_nee_cache(0, 0, 0);      // The last potsition p that can issue a nee. if it’s on an index-matching surface, it can’t issue nee
    bool never_scattered = true;

    while (true) {
        bool scatter = false;
        Spectrum transmittance(1, 1, 1);
        Spectrum pdf_trans(1, 1, 1);

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = infinity<Real>();
        if (vertex_)
            t_hit = distance(vertex_->position, ray.org);
        Real t = infinity<Real>();


        if (current_medium_id != -1) {
            const Medium& current_medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(current_medium, ray.org);
            Spectrum sigma_s = get_sigma_s(current_medium, ray.org);
            Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);

            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            t = log(1 - u) / (-sigma_t[0]);

            // if scatter in the medium
            if (t < t_hit) {
                scatter = true;
				never_scattered = false;
                // pdf eqn (11)
                pdf_trans = sigma_t * exp(-sigma_t * t);
                transmittance = exp(-sigma_t * t);     // transmittance from ray origin to scatter point
                ray.org = ray.org + ray.dir * t;

				pos_nee_cache = ray.org;                  // store this scatter point
            }
            // hit surface 
            else {
                scatter = false;
                pdf_trans = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                // avoid self intersection
                assert(vertex_.has_value());
                ray = Ray{ vertex_->position, ray.dir, get_intersection_epsilon(scene), infinity<Real>() };
            }
        }
        tp = tp * (transmittance / pdf_trans);
		pdf_multi_trans *= pdf_trans[0];

        // nee contribution
        if (scatter) {
            const Medium& current_medium = scene.media[current_medium_id];
            Spectrum sigma_s = get_sigma_s(current_medium, ray.org);
			Spectrum nee_contrib = nee_medium(scene, ray.org, current_medium_id, rng, ray_diff, depth, -ray.dir);
			radiance += tp * sigma_s * nee_contrib;
        }

        // latter part of eqn(14), contribution from hitting emissive geometry surface
        if (!scatter && vertex_) {
			// if hitting the light source.
            if (is_light(scene.shapes[vertex_->shape_id])) {
                Spectrum Le = emission(*vertex_, -ray.dir, scene);
                // This is the only way we can see the light source, so no need MIS
                if (never_scattered) {
                    radiance += tp * Le;
                }
                else {
                    // MIS for hitting light source directly
                    Vector3 light_point = vertex_->position;
                    int light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal point_on_light{ light_point, vertex_->geometric_normal };
                    Real pdf_nee = light_pmf(scene, light_id) * pdf_point_on_light(light, point_on_light, pos_nee_cache, scene);
                    Real G = (abs(dot(-normalize(point_on_light.position - pos_nee_cache), point_on_light.normal))) / (distance_squared(pos_nee_cache, point_on_light.position));
                    Real pdf_phase = pdf_dir_latest * pdf_multi_trans * G;
                    Real w = (pdf_phase * pdf_phase) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

                    radiance += w * tp * Le;
                }
            }
        }

        // reaching max depth
        if (depth == scene.options.max_depth - 1 && scene.options.max_depth != -1)
            break;

        // if hit surface
        if (!scatter && vertex_) {
            // if the surface is not geometric: hitting another participating medium.
            if (vertex_->material_id == -1) {
                // index matching interface, skip through it
                current_medium_id = update_medium_id(*vertex_, ray, current_medium_id);
                // 2/8/2026: THIS IS SOOOOOOOOO  IMPORTANT FOR VACUUM entering participating media.
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                depth++;
                continue;
            }
        }

        // sample next dir & update path tp
        if (scatter){
            PhaseFunction pf = get_phase_function(scene.media[current_medium_id]);
            Vector2 pf_rnd_param{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            std::optional<Vector3> nxt_dir = sample_phase_function(pf, -ray.dir, pf_rnd_param);
            assert(nxt_dir.has_value());
            Real dir_pdf = pdf_sample_phase(pf, -ray.dir, *nxt_dir);
			pdf_dir_latest = dir_pdf;
			pdf_multi_trans = 1;        // reset

            Spectrum rho = eval(pf, *nxt_dir, -ray.dir);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            tp = tp * (rho * sigma_s / dir_pdf);
            ray.dir = *nxt_dir;
        }
        else {  // hit surface, do nothing for now
            break;
        }

        // russian roulette
        Real rr_prob = 1;
        if (depth >= scene.options.rr_depth) {
            rr_prob = min(tp[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
                break;
            else
                tp = tp / rr_prob;
        }
        depth++;
    }
    return radiance;
}


// assume point is on a bsdf material with mat as material
// wi is pointing to camera path
Spectrum nee_bsdf(const Scene& scene, const Vector3& point, int current_medium_id, pcg32_state& rng, RayDifferential& ray_diff, 
            int bounce, const Vector3& wi, const Material& mat, const PathVertex& pathvert, const TexturePool& texturePool) {

    Vector3 p = point;  // copy for later update while walking the ray
    // sample light
    Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light& light = scene.lights[light_id];
    PointAndNormal p_prime =
        sample_point_on_light(light, p, light_uv, shape_w, scene);

    // Compute transmittance to light. Skip through index-matching shapes.
    Spectrum T_light(1, 1, 1);
    int shadow_bounce = 0;
    Real pdf_trans_dir = 1;     // for MIS
    int shadow_medium_id = current_medium_id;
    Real pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, p, scene);

    // this loop is only to update the transmittance && pdf of the path from shading point to the sampled light point.
    while (1) {
        Ray shadow_ray = Ray{ p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime.position, p) };

        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }

        // account for transmittance to next_t
        if (shadow_medium_id != -1) {
            const Medium& shadow_medium = scene.media[shadow_medium_id];
            Spectrum sigma_a = get_sigma_a(shadow_medium, shadow_ray.org);
            Spectrum sigma_s = get_sigma_s(shadow_medium, shadow_ray.org);
            Spectrum sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            pdf_trans_dir *= exp(-sigma_t[0] * next_t);
        }

        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
            // something is blocking, is it opaque surface ?
            if (isect->material_id != -1) {
                // hit opaque surface, blocked
                return make_zero_spectrum();
            }

            // hit index-matching surface, we want to pass through -- this introduces one extra connection vertex
            shadow_bounce++;
            if (scene.options.max_depth != -1 && bounce + shadow_bounce + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }

            shadow_medium_id = update_medium_id(*isect, shadow_ray, shadow_medium_id);
            p = p + next_t * shadow_ray.dir;
        }
    }
    if (T_light[0] > 0) {
        Real G = (abs(dot(-normalize(p_prime.position - point), p_prime.normal))) / (distance_squared(point, p_prime.position));
        Vector3 dir_light = normalize(p_prime.position - point);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        Spectrum f = eval(mat, wi, dir_light, pathvert, texturePool);   // include projection cos.
        Real pdf_bsdf = pdf_sample_bsdf(mat, wi, dir_light, pathvert, texturePool);
        pdf_bsdf *= G;  // solid angle to area space

        Spectrum contrib = T_light * f * G * L / pdf_light;
        Real w = (pdf_light * pdf_light) / (pdf_light * pdf_light + pdf_bsdf * pdf_bsdf);
        return w * contrib;
    }
    return make_zero_spectrum();
}



// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    // MIS is only for direct illumination
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);


    int current_medium_id = scene.camera.medium_id;
    Spectrum tp = Vector3{ 1, 1, 1 };   // path throughput
    Spectrum radiance = make_zero_spectrum();
    int depth = 0;

    /*
        For pdf_nee and pdf_phase, when caculating the MIS weight, the only 2 things I care about are,
        the probability of sampling this light point (getting here) from the last scatter point using phase function sampling:
            I need probability of sampling this dir from phase function, and prob. of sampling the ray length t.
        And
        the probability of sampling this light_point (getting here) from the last scatter point using nee:
            I need probability of sampling this light point, which is pdf_light
    */
    Real pdf_dir_latest = 0;           // In solid angle measure: the pdf of the latest direction sampling (phase or bsdf)
    Real pdf_multi_trans = 1;          // The product PDF of transmittance sampling going through several index-matching surfaces
                                       // from the last phase function sampling or bsdf sampling
    Vector3 pos_nee_cache(0, 0, 0);    // nee_p_cache: The last potsition p that can issue a nee. if it’s on an index-matching surface, it can’t issue nee
    bool never_scattered = true;       // including volume scattering or hitting bsdf geometric surface and "scatter"

    while (true) {
        Vector3 wi = -ray.dir;
        bool scatter = false;
        Spectrum transmittance(1, 1, 1);
        Spectrum pdf_trans(1, 1, 1);

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = infinity<Real>();
        if (vertex_)
            t_hit = distance(vertex_->position, ray.org);
        Real t = infinity<Real>();


        if (current_medium_id != -1) {
            const Medium& current_medium = scene.media[current_medium_id];
            Spectrum sigma_a = get_sigma_a(current_medium, ray.org);
            Spectrum sigma_s = get_sigma_s(current_medium, ray.org);
            Spectrum sigma_t = sigma_a + sigma_s;
            Real u = next_pcg32_real<Real>(rng);

            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            t = log(1 - u) / (-sigma_t[0]);

            // if scatter in the medium
            if (t < t_hit) {
                scatter = true;
                never_scattered = false;
                // pdf eqn (11)
                pdf_trans = sigma_t * exp(-sigma_t * t);
                transmittance = exp(-sigma_t * t);     // transmittance from ray origin to scatter point
                ray.org = ray.org + ray.dir * t;

                pos_nee_cache = ray.org;                  // store this scatter point
            }
            // hit surface 
            else {
                scatter = false;
                pdf_trans = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                // avoid self intersection
                assert(vertex_.has_value());
                ray = Ray{ vertex_->position, ray.dir, get_intersection_epsilon(scene), infinity<Real>() };
            }
        }
        tp = tp * (transmittance / pdf_trans);
        pdf_multi_trans *= pdf_trans[0];

        // volume scattering nee contribution
        if (scatter) {
            const Medium& current_medium = scene.media[current_medium_id];
            Spectrum sigma_s = get_sigma_s(current_medium, ray.org);
            Spectrum nee_contrib = nee_medium(scene, ray.org, current_medium_id, rng, ray_diff, depth, wi);
            radiance += tp * sigma_s * nee_contrib;
        }

        // latter part of eqn(14), contribution from hitting emissive geometry surface
        if (!scatter && vertex_) {
            // if hitting the light source.
            if (is_light(scene.shapes[vertex_->shape_id])) {
                Spectrum Le = emission(*vertex_, wi, scene);
                // This is the only way we can see the light source, so no need MIS
                if (never_scattered) {
                    radiance += tp * Le;
                }
                else {
                    // MIS for hitting light source directly
                    Vector3 light_point = vertex_->position;
                    int light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal point_on_light{ light_point, vertex_->geometric_normal };
                    Real pdf_nee = pdf_point_on_light(light, point_on_light, pos_nee_cache, scene);
                    Real G = (abs(dot(-normalize(point_on_light.position - pos_nee_cache), point_on_light.normal))) / (distance_squared(pos_nee_cache, point_on_light.position));
                    Real pdf_phase = pdf_dir_latest * pdf_multi_trans * G;
                    Real w = (pdf_phase * pdf_phase) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

                    radiance += w * tp * Le;
                }
            }
        }

        // reaching max depth
        if (depth == scene.options.max_depth - 1 && scene.options.max_depth != -1)
            break;

        // if hit surface
        if (!scatter && vertex_) {
            // if the surface is not geometric: hitting another participating medium.
            if (vertex_->material_id == -1) {
                // index matching interface, skip through it
                current_medium_id = update_medium_id(*vertex_, ray, current_medium_id);
                // 2/8/2026: THIS IS SOOOOOOOOO  IMPORTANT FOR VACUUM entering participating media.
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                depth++;
                continue;
            }
        }

        // sample next dir & update path tp
        if (scatter) {
            PhaseFunction pf = get_phase_function(scene.media[current_medium_id]);
            Vector2 pf_rnd_param{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            std::optional<Vector3> nxt_dir = sample_phase_function(pf, wi, pf_rnd_param);
            assert(nxt_dir.has_value());
            Real dir_pdf = pdf_sample_phase(pf, wi, *nxt_dir);
            pdf_dir_latest = dir_pdf;
            pdf_multi_trans = 1;        // reset

            Spectrum rho = eval(pf, *nxt_dir, wi);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            tp = tp * (rho * sigma_s / dir_pdf);
            ray.dir = *nxt_dir;
        }
		// hit geometric surface:
		// 1. calculate dir lighting contrib on geometric surface
        // 2. bsdf sampling for next dir
        else if(!scatter && vertex_){  
            Vector3 geo_hit_point = vertex_->position;
            // nee
            Material mat = scene.materials[vertex_->material_id];
            Vector3 contrib_nee = nee_bsdf(scene, geo_hit_point, current_medium_id, rng, ray_diff, depth, wi, mat, *vertex_, scene.texture_pool);
            radiance += tp * contrib_nee;

            // bsdf sampling
            Vector2 bsdf_rnd_param_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                    wi,
                    *vertex_,
                    scene.texture_pool,
                    bsdf_rnd_param_uv,
                    bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                break;
            }
            Vector3 sampled_wo = bsdf_sample_.value().dir_out;
            Real pdf_bsdf = pdf_sample_bsdf(mat, wi, sampled_wo, *vertex_, scene.texture_pool);
            Spectrum f = eval(mat, wi, sampled_wo, *vertex_, scene.texture_pool);
            tp = tp * f / pdf_bsdf;

            // update some var related to volume medium calculation
            // hitting a geometric surface can be regarded as "scattering happens"
            never_scattered = false;
            pdf_multi_trans = 1;        // reset
            /*
            *  here we sampled direction, can be used for any later potential MIS:
                    Real pdf_phase = dir_pdf_latest * multi_trans_pdf * G;
                    Real w = (pdf_phase * pdf_phase) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

                become ->
                   Real pdf_phase = pdf_bsdf * G;       // actually pdf_bsdf now
                   Real w = (pdf_phase * pdf_phase) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            */
            pdf_dir_latest = pdf_bsdf;  // here we sampled direction, can be used for any later potential MIS
            // offset hit point to avoid shadow acne
            Vector3 N = dot(-ray.dir, (*vertex_).geometric_normal) > 0 ? (*vertex_).geometric_normal : -(*vertex_).geometric_normal;
            pos_nee_cache = (*vertex_).position + N * get_intersection_epsilon(scene);

            // update ray info
            ray = Ray{(*vertex_).position, sampled_wo, get_intersection_epsilon(scene), infinity<Real>() };
            current_medium_id = update_medium_id((*vertex_), ray, current_medium_id);
                                       
        }

        // russian roulette
        Real rr_prob = 1;
        if (depth >= scene.options.rr_depth) {
            rr_prob = min(tp[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
                break;
            else
                tp = tp / rr_prob;
        }
        depth++;
    }
    return radiance;
}

// next event estimation for participating medium, it is used for heterogeneous medium
// wi: output ray, pointing towards the camera path 
Spectrum nee_medium_heter(const Scene& scene, const Vector3& point, int current_medium_id, pcg32_state& rng, RayDifferential& ray_diff, int bounce, const Vector3& wi) {
    Vector3 p = point;  // copy for later update while walking the ray
    // sample light
    Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light& light = scene.lights[light_id];
    PointAndNormal p_prime =
        sample_point_on_light(light, point, light_uv, shape_w, scene);

    // Compute transmittance to light. Skip through index-matching shapes.
    Spectrum T_light(1, 1, 1);
    int shadow_bounce = 0;
    Spectrum pdf_trans_nee = make_const_spectrum(1);     // for MIS
    Spectrum pdf_trans_dir = make_const_spectrum(1);     // for MIS
    int shadow_medium_id = current_medium_id;
    Real pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, point, scene);

    // this loop is only to update the transmittance && pdf of the path from shading point to the sampled light point.
    while (1) {
        Ray shadow_ray = Ray{ p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime.position, p) };

        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }

        // account for transmittance to next_t
        if (shadow_medium_id != -1) {
            const Medium& med = scene.media[shadow_medium_id];
            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Spectrum sigma_m = get_majorant(med, shadow_ray);       // majorant
            Real accum_t = 0;
            int iteration = 0;

            while (1) {
                if (sigma_m[channel] <= 0)
                    break;
                if (iteration >= scene.options.max_null_collisions) 
                    break;
                
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_m[channel];
                Real dt = next_t - accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    // didn't hit surface, so this is  a null-scattering event
                    Spectrum sigma_t = get_sigma_a(med, p) + get_sigma_s(med, p);
                    T_light *= exp(-sigma_m * t) * (sigma_m - sigma_t) / max(sigma_m);
                    pdf_trans_nee *= exp(-sigma_m * t) * sigma_m / max(sigma_m);
                    Spectrum real_prob = sigma_t / sigma_m;
                    pdf_trans_dir *= exp(-sigma_m * t) * sigma_m * (1 - real_prob) / max(sigma_m);
                    p = p + t * shadow_ray.dir;
                    if (max(T_light) <= 0) {    // optimization for places where sigma_n = 0
                        break;
                    }
                }
                else {
                    // hit the surface
                    T_light *= exp(-sigma_m * dt);
                    pdf_trans_dir *= exp(-sigma_m * dt);
                    pdf_trans_nee *= exp(-sigma_m * dt);
                    p = p + dt * shadow_ray.dir;
                    break;
                }
                iteration++;
            }
            
        }

        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
            // something is blocking, is it opaque surface ?
            if (isect->material_id != -1) {
                // hit opaque surface, blocked
                return make_zero_spectrum();
            }

            // hit index-matching surface, we want to pass through -- this introduces one extra connection vertex
            shadow_bounce++;
            if (scene.options.max_depth != -1 && bounce + shadow_bounce + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }

            shadow_medium_id = update_medium_id(*isect, shadow_ray, shadow_medium_id);
            Vector3 N = dot(shadow_ray.dir, isect->geometric_normal) > 0 ?
                isect->geometric_normal : -isect->geometric_normal;
            p = isect->position + N * get_intersection_epsilon(scene);
        }
    }
    if (T_light[0] > 0 || T_light[1] > 0 || T_light[2] > 0) {
        Real G = (abs(dot(-normalize(p_prime.position - point), p_prime.normal))) / (distance_squared(point, p_prime.position));
        Vector3 dir_light = normalize(p_prime.position - point);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        pdf_light *= avg(pdf_trans_nee);

        PhaseFunction pf = get_phase_function(scene.media[current_medium_id]);
        Spectrum rho = eval(pf, wi, dir_light);
        Spectrum contrib = T_light * rho * G * L / pdf_light;
        // MIS: it’s also possible
        // that a phase function sampling + multiple exponential sampling  will reach the light source.
        // We also need to multiply with G to convert phase function PDF to area measure.
        Real pdf_phase = pdf_sample_phase(pf, wi, dir_light) * G * avg(pdf_trans_dir);

        Real w = (pdf_light * pdf_light) / (pdf_light * pdf_light + pdf_phase * pdf_phase);
        return w * contrib;
    }

    return make_zero_spectrum();
}

// assume point is on a bsdf material with mat as material, this is used for heterogeneous volpath
// wi is pointing to camera path
Spectrum nee_bsdf_heter(const Scene& scene, const Vector3& point, int current_medium_id, pcg32_state& rng, RayDifferential& ray_diff,
    int bounce, const Vector3& wi, const Material& mat, const PathVertex& pathvert, const TexturePool& texturePool) {

    Vector3 p = point;  // copy for later update while walking the ray
    // sample light
    Vector2 light_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
    Real light_w = next_pcg32_real<Real>(rng);
    Real shape_w = next_pcg32_real<Real>(rng);
    int light_id = sample_light(scene, light_w);
    const Light& light = scene.lights[light_id];
    PointAndNormal p_prime =
        sample_point_on_light(light, p, light_uv, shape_w, scene);

    // Compute transmittance to light. Skip through index-matching shapes.
    Spectrum T_light(1, 1, 1);
    int shadow_bounce = 0;
    Spectrum pdf_trans_dir = make_const_spectrum(1);
    Spectrum pdf_trans_nee = make_const_spectrum(1);
    int shadow_medium_id = current_medium_id;
    Real pdf_light = light_pmf(scene, light_id) * pdf_point_on_light(light, p_prime, point, scene);

    // this loop is only to update the transmittance && pdf of the path from shading point to the sampled light point.
    while (1) {
        Ray shadow_ray = Ray{ p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) *
                                    distance(p_prime.position, p) };

        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }

        // account for transmittance to next_t
        if (shadow_medium_id != -1) {
            const Medium &med = scene.media[shadow_medium_id];
            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Spectrum sigma_m = get_majorant(med, shadow_ray);       // majorant
            Real accum_t = 0;
            int iteration = 0;

            while (1) {
                if (sigma_m[channel] <= 0)
                    break;
                if (iteration >= scene.options.max_null_collisions)
                    break;

                Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_m[channel];
                Real dt = next_t - accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    // didn't hit surface, so this is  a null-scattering event
                    Spectrum sigma_t = get_sigma_a(med, p) + get_sigma_s(med, p);
                    T_light *= exp(-sigma_m * t) * (sigma_m - sigma_t) / max(sigma_m);
                    pdf_trans_nee *= exp(-sigma_m * t) * sigma_m / max(sigma_m);
                    Spectrum real_prob = sigma_t / sigma_m;
                    pdf_trans_dir *= exp(-sigma_m * t) * sigma_m * (1 - real_prob) / max(sigma_m);
                    p = p + t * shadow_ray.dir;
                    if (max(T_light) <= 0) {    // optimization for places where sigma_n = 0
                        break;
                    }
                }
                else {
                    // hit the surface
                    T_light *= exp(-sigma_m * dt);
                    pdf_trans_dir *= exp(-sigma_m * dt);
                    pdf_trans_nee *= exp(-sigma_m * dt);
                    p = p + dt * shadow_ray.dir;
                    break;
                }
                iteration++;
            }
            
        }
        if (!isect) {
            // Nothing is blocking, we’re done
            break;
        }
        else {
            // something is blocking, is it opaque surface ?
            if (isect->material_id != -1) {
                // hit opaque surface, blocked
                return make_zero_spectrum();
            }

            // hit index-matching surface, we want to pass through -- this introduces one extra connection vertex
            shadow_bounce++;
            if (scene.options.max_depth != -1 && bounce + shadow_bounce + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }


            shadow_medium_id = update_medium_id(*isect, shadow_ray, shadow_medium_id);

            Vector3 N = dot(shadow_ray.dir, isect->geometric_normal) > 0 ?
                isect->geometric_normal : -isect->geometric_normal;
            p = isect->position + N * get_intersection_epsilon(scene);
           //  p = p + next_t * shadow_ray.dir;
        }
    }
    if (T_light[0] > 0 || T_light[1] > 0 || T_light[2] > 0) {
        Real G = (abs(dot(-normalize(p_prime.position - point), p_prime.normal))) / (distance_squared(point, p_prime.position));
        Vector3 dir_light = normalize(p_prime.position - point);
        Spectrum L = emission(light, -dir_light, Real(0), p_prime, scene);
        Spectrum f = eval(mat, wi, dir_light, pathvert, texturePool);   // include projection cos.
        pdf_light *= avg(pdf_trans_nee);

        Real pdf_bsdf = pdf_sample_bsdf(mat, wi, dir_light, pathvert, texturePool);
        pdf_bsdf *= G;  // solid angle to area space
        pdf_bsdf *= avg(pdf_trans_dir);
        Spectrum contrib = T_light * f * G * L / pdf_light;
           
        Real w = (pdf_light * pdf_light) / (pdf_light * pdf_light + pdf_bsdf * pdf_bsdf);
        return w * contrib;
    }
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);

    int current_medium_id = scene.camera.medium_id;
    Spectrum tp = Vector3{ 1, 1, 1 };   // path throughput
    Spectrum radiance = make_zero_spectrum();
    int depth = 0;

    Spectrum  pdf_dir_latest(0, 0, 0);          // In solid angle measure: the pdf of the latest direction sampling (phase or bsdf)
    Spectrum  pdf_multi_trans(1,1,1);           // The product PDF of transmittance sampling going through several index-matching surfaces
                                                // from the last phase function sampling or bsdf sampling
    Vector3 pos_nee_cache(0, 0, 0);             // The last potsition p that can issue a nee. if it’s on an index-matching surface, it can’t issue nee

    // hw2 Slides page 20:
    // To combine phase function sampling with the next event estimation, we
    // need the PDFs for the ratio tracking (the sequence of decisions we make to evaluate either the recursive term
    // or the source term) during free - flight sampling, and the PDFs of free - flight sampling(the real / fake particle
    // events) during next event estimation.
    Spectrum pdf_multi_nee = make_const_spectrum(Real(1)); 

    bool never_scattered = true;                // Including volume scattering or hitting bsdf geometric surface and "scatter"


    while (1) {
        bool scatter = false;
        Vector3 wi = -ray.dir;

        Spectrum transmittance = make_const_spectrum(Real(1));
        Spectrum pdf_trans_dir = make_const_spectrum(Real(1));  // PDF for free-flight sampling
        Spectrum pdf_trans_nee = make_const_spectrum(Real(1));  // PDF for next event estimation

        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Real t_hit = infinity<Real>();
        if (vertex_)
            t_hit = distance(vertex_->position, ray.org);
        Real t = infinity<Real>();


        if (current_medium_id != -1) {
            const Medium& med = scene.media[current_medium_id];
            Spectrum sigma_m = get_majorant(med, ray);  // majorant

            // Sample a channel for sampling
            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;

            // homogenized free-flight sampling
            while (1) {
                if (sigma_m[channel] <= 0) 
                    break;

                if (iteration >= scene.options.max_null_collisions)
                    break;

                t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_m[channel];
                Real dt = t_hit - accum_t;
                // Update accumulated distance
                accum_t = min(accum_t + t, t_hit);
                if (t < dt) {   // haven’t reached the surface
                    // sample from real / fake particle events
                    Vector3 pos_current_p = ray.org + accum_t * ray.dir;
                    Spectrum sigma_t = get_sigma_a(med, pos_current_p) + get_sigma_s(med, pos_current_p);
                    Spectrum real_prob = sigma_t / sigma_m;

                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        //  hit a "real" particle
                        scatter = true;
                        transmittance *= exp(-sigma_m * t) / max(sigma_m);
                        pdf_trans_dir *= exp(-sigma_m * t) * sigma_m * real_prob / max(sigma_m);
                        // don't need to account for pdf_trans_nee since we scatter
                        ray.org = ray.org + accum_t * ray.dir;
                        break;
                    }
                    else {
                        // hit a "fake" particle
                        transmittance *= exp(-sigma_m * t) * (sigma_m - sigma_t) / max(sigma_m);
                        pdf_trans_dir *= exp(-sigma_m * t) * sigma_m * (1 - real_prob) / max(sigma_m);
                        pdf_trans_nee *= exp(-sigma_m * t) * sigma_m / max(sigma_m);
                    }
                }
                // reach the surface
                else {
                    transmittance *= exp(-sigma_m * dt);
                    pdf_trans_dir *= exp(-sigma_m * dt);
                    pdf_trans_nee *= exp(-sigma_m * dt);
                    
                    // important
                    ray.org = ray.org + t_hit * ray.dir;
                    break;
                }
                iteration++;
            }
        }
        if (current_medium_id == -1 && !scatter && !vertex_) {
            break;   // terminate path
        }

        tp *= transmittance / avg(pdf_trans_dir);
        pdf_multi_trans *= pdf_trans_dir;
        pdf_multi_nee *= pdf_trans_nee;

        // volume scattering nee contribution
        if (scatter) {
            never_scattered = false;
            const Medium& current_medium = scene.media[current_medium_id];
            Spectrum sigma_s = get_sigma_s(current_medium, ray.org);
            Spectrum nee_contrib = nee_medium_heter(scene, ray.org, current_medium_id, rng, ray_diff, depth, wi);
            radiance += tp * sigma_s * nee_contrib;
        }

        // latter part of eqn(14), contribution from hitting emissive geometry surface
        if (!scatter && vertex_) {
            // if hitting the light source.
            if (is_light(scene.shapes[vertex_->shape_id])) {
                Spectrum Le = emission(*vertex_, wi, scene);
                // This is the only way we can see the light source, so no need MIS
                if (never_scattered) {
                    radiance += tp * Le;
                }
                else {
                    // MIS for hitting light source directly
                    Vector3 light_point = vertex_->position;
                    int light_id = get_area_light_id(scene.shapes[(*vertex_).shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal point_on_light{ light_point, vertex_->geometric_normal };

                    Spectrum pdf_nee = pdf_multi_nee * pdf_point_on_light(light, point_on_light, pos_nee_cache, scene);

                    Real G = (abs(dot(-normalize(point_on_light.position - pos_nee_cache), point_on_light.normal))) / (distance_squared(pos_nee_cache, point_on_light.position));
                    Spectrum pdf_phase = pdf_dir_latest * pdf_multi_trans * G;
                    Spectrum w = (pdf_phase * pdf_phase) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

                    radiance += w * tp * Le;
                }
            }
        }

        // reaching max depth
        if (depth == scene.options.max_depth - 1 && scene.options.max_depth != -1)
            break;

        // if hit surface
        if (!scatter && vertex_) {
            // if the surface is not geometric: hitting another participating medium.
            if (vertex_->material_id == -1) {
                // index matching interface, skip through it
                current_medium_id = update_medium_id(*vertex_, ray, current_medium_id);
                // 2/8/2026: THIS IS SOOOOOOOOO  IMPORTANT FOR VACUUM entering participating media.
                ray.org = vertex_->position + ray.dir * get_intersection_epsilon(scene);
                depth++;
                continue;
            }
        }

        // sample next dir & update path tp
        if (scatter) {
            PhaseFunction pf = get_phase_function(scene.media[current_medium_id]);
            Vector2 pf_rnd_param{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            std::optional<Vector3> nxt_dir = sample_phase_function(pf, wi, pf_rnd_param);
            assert(nxt_dir.has_value());
            Real dir_pdf = pdf_sample_phase(pf, wi, *nxt_dir);
            pdf_dir_latest = Vector3(dir_pdf, dir_pdf, dir_pdf);        // 2/10/2026: ?
            pdf_multi_trans = make_const_spectrum(1);        // reset
            pdf_multi_nee = make_const_spectrum(1);          // reset ?
            pos_nee_cache = ray.org;

            Spectrum rho = eval(pf, *nxt_dir, wi);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium_id], ray.org);
            tp = tp * (rho * sigma_s / dir_pdf);
            ray.dir = *nxt_dir;
        }
        else if (!scatter && vertex_) {
            Vector3 geo_hit_point = vertex_->position;
            // nee
            const Material &mat = scene.materials[vertex_->material_id];

            Vector3 contrib_nee = nee_bsdf_heter(scene, geo_hit_point, current_medium_id, rng, ray_diff, depth, wi, mat, *vertex_, scene.texture_pool);

            radiance += tp * contrib_nee;

            // bsdf sampling
            Vector2 bsdf_rnd_param_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
            Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
            std::optional<BSDFSampleRecord> bsdf_sample_ =
                sample_bsdf(mat,
                    wi,
                    *vertex_,
                    scene.texture_pool,
                    bsdf_rnd_param_uv,
                    bsdf_rnd_param_w);
            if (!bsdf_sample_) {
                break;
            }
            Vector3 sampled_wo = bsdf_sample_.value().dir_out;
            Real pdf_bsdf = pdf_sample_bsdf(mat, wi, sampled_wo, *vertex_, scene.texture_pool);
            Spectrum f = eval(mat, wi, sampled_wo, *vertex_, scene.texture_pool);
            tp = tp * f / pdf_bsdf;

            // update some var related to volume medium calculation
            // hitting a geometric surface can be regarded as "scattering happens"
            never_scattered = false;
            pdf_multi_trans = make_const_spectrum(1);        // reset
            pdf_multi_nee = make_const_spectrum(1);          // reset

            pdf_dir_latest = make_const_spectrum(pdf_bsdf);  // here we sampled direction, can be used for any later potential MIS
            // offset hit point to avoid shadow acne
            Vector3 N = dot(-ray.dir, (*vertex_).geometric_normal) > 0 ? (*vertex_).geometric_normal : -(*vertex_).geometric_normal;
            pos_nee_cache = (*vertex_).position + N * get_intersection_epsilon(scene);

            // update ray info
            ray = Ray{ (*vertex_).position, sampled_wo, get_intersection_epsilon(scene), infinity<Real>() };
            current_medium_id = update_medium_id((*vertex_), ray, current_medium_id);
        }

        // russian roulette
        Real rr_prob = 1;
        if (depth >= scene.options.rr_depth) {
            rr_prob = min(tp[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob)
                break;
            else
                tp = tp / rr_prob;
        }
        depth++;
    
    }
    return radiance;
}
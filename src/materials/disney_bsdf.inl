#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool isInside = dot(vertex.geometric_normal, dir_in) <= 0;

    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Spectrum f_diffuse = eval_op::operator()(DisneyDiffuse{
        bsdf.base_color,
        bsdf.roughness,
        bsdf.subsurface });

    Spectrum f_sheen = eval_op::operator()(DisneySheen{
        bsdf.base_color,
		bsdf.sheen_tint });

	Spectrum f_clearcoat = eval_op::operator()(DisneyClearcoat{
        bsdf.clearcoat_gloss });

    Spectrum f_glass = eval_op::operator()(DisneyGlass{
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic,
        bsdf.eta });


    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);

    // weights 
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = 1 - specular_transmission * (1 - metallic);
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real w_sheen = (1 - metallic) * sheen;

	Vector3 h = normalize(dir_in + dir_out);
    if(!reflect) {
        h = normalize(dir_in + dir_out * eta);
	}

	Spectrum f_metal = Spectrum(0, 0, 0);
    // ******** f_metal begins
    {
        Real R0 = (eta - 1) * (eta - 1) / ((eta + 1) * (eta + 1));
        Real lum = luminance(base_color);
        Spectrum Ctint = lum > 0 ? base_color / lum : make_const_spectrum(1);
        Spectrum Ks = (1 - specular_tint) + specular_tint * Ctint;
	    Spectrum C0 = specular * R0 * (1 - metallic) * Ks + metallic * base_color;
	    Spectrum fresnel_hat = C0 + (1 - C0) * pow(1 - abs(dot(h, dir_out)), 5);

        Vector3 w_in = dir_in;
        Vector3 w_out = dir_out;

        Spectrum F_m = fresnel_hat;

        Real D_m = D_metal(h, ax, ay, frame);
        Real G_m = G(w_in, ax, ay, frame) * G(w_out, ax, ay, frame);

        f_metal = F_m * G_m * D_m / (4 * abs(dot(frame.n, w_in)));

        if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
            // No light below the surface
            f_metal = Spectrum(0, 0, 0);
        }
    }
    // ******** f_metal ends

    if (isInside) {
		f_clearcoat = Spectrum(0, 0, 0);
		f_diffuse = Spectrum(0, 0, 0);
		f_metal = Spectrum(0, 0, 0);
		f_sheen = Spectrum(0, 0, 0);
    }

    // final f
	Spectrum f_disney = 
        w_diffuse * f_diffuse +
        w_sheen * f_sheen +
        w_metal * f_metal +
        w_clearcoat * f_clearcoat +
        w_glass * f_glass;

	return f_disney;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    /*
    do not obtain pdf like below, hardcode instead

    Real pdf_glass = pdf_sample_bsdf_op::operator()(DisneyGlass{
        bsdf.base_color,
        bsdf.roughness,
        bsdf.anisotropic,
		bsdf.eta });

    ...

    */

	bool isInside = dot(vertex.geometric_normal, dir_in) < 0;

    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // parameters
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);

    Vector3 w_in = dir_in;
    Vector3 w_out = dir_out;

	Vector3 h = normalize(dir_in + dir_out);
    if (!reflect) {
        h = normalize(dir_in + dir_out * eta);
    }

	// invert half-vector if it's on  the other side
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    Real n_dot_in = dot(frame.n, dir_in);

	// Glass PDF
    Real pdf_glass = 0;
    Real h_dot_in = dot(h, dir_in);

    Real Rs = (dot(h, w_in) - eta * dot(h, w_out)) / (dot(h, w_in) + eta * dot(h, w_out));
    Real Rp = (eta * dot(h, w_in) - dot(h, w_out)) / (eta * dot(h, w_in) + dot(h, w_out));
    Real F_g = fresnel_dielectric(h_dot_in, eta);

    Real D_g = D_metal(h, ax, ay, frame);
    Real G_g = G(w_in, ax, ay, frame);
    if (reflect) {
        pdf_glass = (F_g * D_g * G_g) / (4 * fabs(dot(frame.n, dir_in)));
    }
    else {
        Real h_dot_out = dot(h, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        pdf_glass = (1 - F_g) * D_g * G_g * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }

    // Metal PDF
    Real D_m = D_metal(h, ax, ay, frame);
    Real G_in = G(w_in, ax, ay, frame);
    Real pdf_metal = (G_in * D_m) / (4 * abs(n_dot_in));

    // Diffuse PDF
    Real pdf_diffuse = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;


	// Clearcoat PDF
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
    Vector3 h_local = to_local(frame, h);

    Real Dc = (alpha_g * alpha_g - 1) /
        (c_PI * log(alpha_g * alpha_g) * (1 + (alpha_g * alpha_g - 1) * h_local.z * h_local.z));

    Real pdf_clearcoat = Dc * abs(dot(frame.n, h) / (4 * abs(dot(h, dir_out))));


    // weights 
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = 1 - specular_transmission * (1 - metallic);
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real w_sum = w_diffuse + w_metal + w_clearcoat + w_glass;
    if(w_sum <= 0) 
        return 0;
	
    // normalization
    w_diffuse /= w_sum;
    w_metal /= w_sum;
    w_clearcoat /= w_sum;
    w_glass /= w_sum;
    
    if (isInside) {
		w_diffuse = 0;
		w_metal = 0;
		w_clearcoat = 0;
		w_glass = 1;
    }
    
	// Pdf(w_out | pickedLobe = material) = p_material;
	// Pdf(w_out, pickedLobe = material) = w_material * p_material;
	// pdf(w_out) = sum_over_all_material( w_material' * p_material' )
	Real p = w_diffuse * pdf_diffuse + w_metal * pdf_metal + w_clearcoat * pdf_clearcoat + w_glass * pdf_glass;
    return p;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // parameters
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);


    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);

    // weights 
    Real w_diffuse = (1 - specular_transmission) * (1 - metallic);
    Real w_metal = 1 - specular_transmission * (1 - metallic);
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

	Real w_sum = w_diffuse + w_metal + w_clearcoat + w_glass;
    // normalization
	w_diffuse /= w_sum;
	w_metal /= w_sum;
	w_clearcoat /= w_sum;
	w_glass /= w_sum;

	// if ray is inside the object, only sample glass
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        return sample_bsdf_op::operator()(DisneyGlass{
            bsdf.base_color,
            bsdf.roughness,
            bsdf.anisotropic,
            bsdf.eta });
	}

	// randomly choose a lobe to sample
	Real r = rnd_param_w; 

    if (r < w_diffuse) {
        // sample diffuse
        return sample_bsdf_op::operator()(DisneyDiffuse{
            bsdf.base_color,
            bsdf.roughness,
            bsdf.subsurface });
    }
    else if (r < w_diffuse + w_metal) {
        // sample metal
        return sample_bsdf_op::operator()(DisneyMetal{
           bsdf.base_color,
           bsdf.roughness,
           bsdf.anisotropic });
    }
    else if (r < w_diffuse + w_metal + w_clearcoat) {
        // sample clearcoat
        return sample_bsdf_op::operator()(DisneyClearcoat{
			bsdf.clearcoat_gloss });
    }
    else {
        // sample glass
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, ax, ay, rnd_param_uv);

        Vector3 h = to_world(frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(h, frame.n) < 0) {
            h = -h;
        }
        Real h_dot_in = dot(h, dir_in);

        Vector3 w_in = dir_in;

        // I gave up on using exact fresnel eqn here.
        //Real Rs = (dot(h, w_in) - eta * dot(h, w_out)) / (dot(h, w_in) + eta * dot(h, w_out));
        //Real Rp = (eta * dot(h, w_in) - dot(h, w_out)) / (eta * dot(h, w_in) + dot(h, w_out));
        Real F_g = fresnel_dielectric(h_dot_in, eta);

		// rescale the random number:
		// a = w_diffuse + w_metal + w_clearcoat
        // now r is in [a, 1);
		Real rnd = (r - (w_diffuse + w_metal + w_clearcoat)) / w_glass;

        if (rnd <= F_g) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);
            // set eta to 0 since we are not transmitting
            return BSDFSampleRecord{ reflected, Real(0) /* eta */, roughness };
        }
        else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0) {
                // Total internal reflection
                // This shouldn't really happen, as F will be 1 in this case.
                return {};
            }
            // flip half_vector if needed
            if (h_dot_in < 0) {
                h = -h;
            }
            Real h_dot_out = sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * h;
            return BSDFSampleRecord{ refracted, eta, roughness };
        }
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}

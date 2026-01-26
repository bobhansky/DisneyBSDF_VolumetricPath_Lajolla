#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);

    if(dot(frame.n, dir_in) <= 0 || dot(frame.n, dir_out) <= 0) {
        return make_zero_spectrum();
	};

    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
   
    Vector3 h_local = to_local(frame, h);
    Real R0 = 0.25 / (2.5 * 2.5);
    Real Fc = R0 + (1 - R0) * pow(1 - abs(dot(h, dir_out)), 5);

    Real Dc = (alpha_g * alpha_g - 1) /
        (c_PI * log(alpha_g * alpha_g) * (1 + (alpha_g * alpha_g - 1) * h_local.z * h_local.z));

    auto lambda = [&](const Vector3& w) {
        Vector3 local_w = to_local(frame, w);
        Real res = sqrt(1 + ((local_w.x * 0.25) * (local_w.x * 0.25) + (local_w.y * 0.25) * (local_w.y * 0.25)) / ((local_w.z) * (local_w.z))) - 1;
        res /= 2;
        return res;
        };
    Real Gc = 1.0 / ((1 + lambda(dir_in)) * (1 + lambda(dir_out)));

	Real f_clearcoat = Fc * Gc * Dc / (4 * abs(dot(frame.n, dir_in)));

    return make_const_spectrum(f_clearcoat);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);

    if (dot(frame.n, dir_in) <= 0 || dot(frame.n, dir_out) <= 0) {
        return 0;
    }
    if (dot(frame.n, h) <= 0) {
        return 0;
    }

	Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
	
	Vector3 h_local = to_local(frame, h);

    Real Dc = (alpha_g * alpha_g - 1) /
        (c_PI * log(alpha_g * alpha_g) * (1 + (alpha_g * alpha_g - 1) * h_local.z * h_local.z));

    Real res = Dc * abs(dot(frame.n, h) / (4 * abs(dot(h, dir_out))));

    return res;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

	Real u = rnd_param_uv[0];
	Real v = rnd_param_uv[1];
	Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;

	Real cos_theta = sqrt( (1 - pow(alpha_g * alpha_g, 1.0 - u)) / (1 - alpha_g * alpha_g));
	Real sin_theta = sqrt(1 - cos_theta * cos_theta);
	Real phi = 2 * c_PI * v;

    Vector3 h_local = Vector3{
		sin_theta* cos(phi),
		sin_theta* sin(phi),
		cos_theta
	};

	Vector3 h = to_world(frame, h_local);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, h) * h);

    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, 0, /* roughness */
    };

}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}

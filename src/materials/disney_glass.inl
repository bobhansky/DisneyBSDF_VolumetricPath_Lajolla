#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

	Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
	Vector3 w_in = dir_in;
	Vector3 w_out = dir_out;

	Real D_m = D_metal(h, ax, ay, frame);

	Real Gg = G(w_in, ax, ay, frame) * G(w_out, ax, ay, frame);

    Real h_dot_in = dot(h, w_in);
	Real F_g = fresnel_dielectric(h_dot_in, eta);
    //Real Rs = (dot(h, w_in) - eta * dot(h, w_out)) / (dot(h, w_in) + eta * dot(h, w_out));
    //Real Rp = (eta * dot(h, w_in) - dot(h, w_out)) / (eta * dot(h, w_in) + dot(h, w_out));
	// F_g = (Rs * Rs + Rp * Rp) / 2;

    if (reflect) {
        Spectrum f_glass = base_color * F_g * D_m * Gg / (4 * fabs(dot(frame.n, w_in)));
        return f_glass;
    }
    else {
        Real h_dot_out = dot(h, w_out);
		Spectrum numerator = sqrt(base_color) * (1 - F_g) * D_m * Gg * fabs(h_dot_out * h_dot_in);
		Real denominator = fabs(dot(frame.n, w_in)) * pow(h_dot_in + eta * h_dot_out, 2);
		Spectrum f_glass = numerator / denominator;
		return f_glass;
	}
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);
    Vector3 w_in = dir_in;
    Vector3 w_out = dir_out;


    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(h, dir_in);

	Real F_g = fresnel_dielectric(h_dot_in, eta);
    // Real Rs = (dot(h, w_in) - eta * dot(h, w_out)) / (dot(h, w_in) + eta * dot(h, w_out));
	//Real Rp = (eta * dot(h, w_in) - dot(h, w_out)) / (eta * dot(h, w_in) + dot(h, w_out));
    // F_g = (Rs * Rs + Rp * Rp) / 2;

	Real D_g = D_metal(h, ax, ay, frame);
	Real G_g = G(w_in, ax, ay, frame);
    if (reflect) {
        return (F_g * D_g * G_g) / (4 * fabs(dot(frame.n, dir_in)));
    }
    else {
        Real h_dot_out = dot(h, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F_g) * D_g * G_g * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1 - 0.9 * anisotropic);
    Real ax = max(0.001, roughness * roughness / aspect);
    Real ay = max(0.001, roughness * roughness * aspect);

    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, ax, ay, rnd_param_uv);

    Vector3 h = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(h, dir_in);

    Vector3 w_in = dir_in;

    // I gave up on using exact fresnel eqn here.
    Real h_dot_out = sqrt(1 - (1 - h_dot_in * h_dot_in) / (eta * eta));
	Real F_g = fresnel_dielectric(h_dot_in, eta);
    //Real Rs = (dot(h, w_in) - eta * h_dot_out) / (dot(h, w_in) + eta * h_dot_out);
    //Real Rp = (eta * dot(h, w_in) - h_dot_out) / (eta * dot(h, w_in) + h_dot_out);
    // F_g = (Rs * Rs + Rp * Rp) / 2;


    if (rnd_param_w <= F_g) {
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

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}

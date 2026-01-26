Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // Implementation based on homework1.pdf
	Vector3 h = normalize(dir_in + dir_out);
	Vector3 base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
	Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

	const Vector3 &w_in = dir_in;
    const Vector3 &w_out = dir_out;

	Real F_D90 = 0.5 + 2.0 * roughness * abs(dot(h, w_out)) * abs(dot(h, w_out));
    // functor
    auto F_D = [&](const Vector3 &w) {
        return (1 + (F_D90 - 1) * pow(1 - abs(dot(frame.n, w)), 5));
		};
    
	Spectrum f_baseDiffuse = (base_color / c_PI) * F_D(w_in) * F_D(w_out) * abs(dot(frame.n, w_out));

	Real F_SS90 = roughness * abs(dot(h, w_out)) * abs(dot(h, w_out));
    auto F_SS = [&](const Vector3 &w) {
        return (1 + (F_SS90 - 1) * pow(1 - abs(dot(frame.n, w)), 5));
		};

	Vector3 f_subsurface = 1.25 * (base_color / c_PI) * 
                          (F_SS(w_in) * F_SS(w_out) * 
                          (1 / abs((dot(frame.n, w_in)) + abs(dot(frame.n, w_out)) ) - 0.5) +
							  0.5) * abs(dot(frame.n, w_out));

	Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
	Spectrum f_diffuse = (1 - subsurface) * f_baseDiffuse + subsurface * f_subsurface;
    return f_diffuse;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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

	// retrieve eta and roughness from the DisneyDiffuse BSDF
    Real eta(0);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness /* roughness */ };
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}

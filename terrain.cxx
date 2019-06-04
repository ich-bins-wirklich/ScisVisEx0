// This source code is property of the Computer Graphics and Visualization chair of the
// TU Dresden. Do not distribute! 
// Copyright (C) CGV TU Dresden - All Rights Reserved

#include <cgv/base/node.h> // this should be first include to avoid warning under VS
#include <cgv/base/register.h>
#include <cgv/data/data_view.h>
#include <cgv/math/ftransform.h>
#include <cgv/media/image/image_reader.h>
#include <cgv/media/mesh/simple_mesh.h>
#include <cgv_reflect_types/media/color.h>
#include <cgv_gl/gl/gl.h>
#include <cgv_gl/sphere_renderer.h>
#include <cgv/render/drawable.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv/render/vertex_buffer.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/event_handler.h>
#include <cgv/gui/key_event.h>
#include <random>
#include <queue>

#include <math.h>

class terrain :
	public cgv::base::node,
	public cgv::render::multi_pass_drawable,
	public cgv::gui::event_handler,
	public cgv::gui::provider
{
private:

protected:
	/**@name hierarchy of right triangles*/
	//@{
	/// texture dimension which yields arrays for height, error, and radius of (N+1)x(N+1)
	size_t N;
	/// number of subdivisions along one coordinate axis in a triangle batch (1 corresponds to individual triangles, assumed to be power of 2)
	unsigned subdivide_count;
	/// depth of tree of right triangles 
	short tree_depth;
	/// height value that corresponds to 1 (255 for 8-bit dem textures and 65535 for 16 bit)
	unsigned max_dem_value;
	/// extent of the terrain in world space
	vec3 extent;
	/// (N+1)x(N+1) height values, where the last row and column are replicated
	std::vector<short> heights;
	/// function to convert index pair to linear index
	int index(short x, short y) const
	{
		return y * (N + 1) + x;
	}
	/// type of triangle node in hierarchy
	struct triangle_node
	{
		// texel location in the range [0,N]
		short x, y;
		// length of triangle along short edge of triangles in even levels and along long edge for triangles in odd levels 
		short base_length;
		// triangle orientation in the range [0,7]
		short omega;
	};
	/// overload for direct index computation on triangle node
	int index(const triangle_node& n) const
	{
		return index(n.x,n.y);
	}
	/// return one of the two root triangles covering the rectangular domain
	triangle_node get_root_triangle(short i) const
	{
		triangle_node root;
		root.base_length = short(N);
		root.x = short(N) / 2;
		root.y = short(N) / 2;
		root.omega = 4 * i;
		return root;
	}
	/// check whether node is leaf
	bool is_leaf(const triangle_node& n) const
	{
		return n.base_length <= (int)subdivide_count;
	}
	/// integer implementation of binary logarithm
	static unsigned log2i(unsigned v)
	{
		unsigned targetlevel = 0;
		while (v >>= 1) 
			++targetlevel;
		return targetlevel;
	}
	/// return the level of a given node
	short level(const triangle_node& n) const
	{
		return 2 * log2i(N / n.base_length) + (n.omega&1);
	}

	static vec2 nav_vertex(short base_length, short omega)
	{
		short x;
		short y;

		if (omega % 2 == 0)
		{
			x = omega == 0 ? 1 : omega == 4 ? -1 : 0;
			y = omega == 2 ? 1 : omega == 6 ? -1 : 0;

			x *= base_length;
			y *= base_length;
		} else
		{
			x = omega == 1 || omega == 7 ? base_length / 2 : -base_length / 2;
			y = omega == 1 || omega == 3 ? base_length / 2 : -base_length / 2;
		}

		return vec2(x, y);


		/*

		*********** this solution produces numerically inaccurate values **********

		float alpha = omega * 3.14159 / 4.0;
		short x = (short) (base_length * cos(alpha));
		short y = (short) (base_length * sin(alpha));

		if (omega % 2 == 1)
		{
			x /= sqrt(2.0);
			y /= sqrt(2.0);
		}

		return vec2(x, y);
		*/
	}
	
	/// return triangle node of left child
	static triangle_node left_child(const triangle_node& n)
	{
		triangle_node res;
		/************************************************************************************
		 tasks 3.1a: Compute the left child of the current node. */

		 /*<your_code_here>*/
		vec2 vertex = nav_vertex(n.base_length / 2, (n.omega + 4) % 8);
		res.x = n.x + vertex[0];
		res.y = n.y + vertex[1];

		res.base_length = n.omega % 2 == 0 ? n.base_length : n.base_length / 2;
		res.omega = (n.omega + 3) % 8;
		/************************************************************************************/
		return res;
	}
	/// return triangle node of right child
	static triangle_node right_child(const triangle_node& n)
	{
		triangle_node res;
		/************************************************************************************
		 tasks 3.1a: Compute the right child of the current node. */

		 /*<your_code_here>*/
		vec2 vertex = nav_vertex(n.base_length / 2, (n.omega + 6) % 8);
		res.x = n.x + vertex[0];
		res.y = n.y + vertex[1];

		res.base_length = n.omega % 2 == 0 ? n.base_length : n.base_length / 2;
		res.omega = (n.omega + 5) % 8;
		/************************************************************************************/
		return res;
	}
	/// check whether neighbor of triangle in diamond top is still inside of domain and in this case set nn to diamond neighbor
	bool has_diamond_neighbor(const triangle_node& n, triangle_node& nn) const
	{
		if (n.x == 0 || n.y == 0 || n.x == N || n.y == N)
			return false;
		nn = n;
		nn.omega = (nn.omega + 4) & 7;
		return true;
	}
	/// return the 3d world point of the diamond point of a triangle
	vec3 world_point(const triangle_node& n) const
	{
		return vec3(
			extent(0)*float(n.x) / N - 0.5f*extent(0),
			extent(1)*float(n.y) / N - 0.5f*extent(1),
			extent(2)*heights[index(n)]/max_dem_value);
	}
	//@}

	/**@name adaptation*/
	//@{
	/// per tex the radius of a sphere in world space 
	std::vector<float> radii, radii_cache;
	/// the height error, scaled such that heights range in [0,1] and the world error is computed by multiplying with extent(2)
	std::vector<float> errors, errors_cache;
	///

private:
	/// per texel flag used to keep track of processed nodes
	std::vector<bool> processed;

protected:
	/// adaptation mode enum
	enum AdaptationMode {
		AM_NONE,
		AM_TREE_DEPTH,
		AM_ISOTROPIC_ERROR,
		AM_ANISOTROPIC_ERROR,
		AM_FIRST = AM_NONE,
		AM_LAST = AM_ISOTROPIC_ERROR
	};
	/// selected adaptation mode
	AdaptationMode adaptation_mode;
	/// adaptation depth for AM_TREE_DEPTH
	short adapted_tree_depth;
	/// maximum adaptation depth during error based adaptation
	short max_tree_depth;
	/// pixel threshold used for adaptation
	float pixel_threshold;
	/// 2*tan_of_half_of_fovy/height_in_pixel
	float view_factor;
	/// world location of eye
	vec3 eye_world;
	/// error of root diamond
	float root_error;
	/// radius of root sphere
	float root_radius;
	/// update slider ranges in gui
	void configure_gui()
	{
		if (find_control(max_tree_depth))
			find_control(max_tree_depth)->set("max", tree_depth);
		if (find_control(adapted_tree_depth))
			find_control(adapted_tree_depth)->set("max", tree_depth);
	}
	/// to be called whenever texture resolution N or subdivide count is altered
	void update_tree_depth()
	{
		tree_depth = 0;
		size_t _N = N / subdivide_count;
		while (_N > 1) {
			++tree_depth;
			_N /= 2;
		}
		tree_depth = 2 * tree_depth - 1;
		max_tree_depth = tree_depth;
		update_member(&max_tree_depth);
		if (adapted_tree_depth > tree_depth) {
			adapted_tree_depth = tree_depth;
			update_member(&adapted_tree_depth);
		}
		configure_gui();
	}
	/// set a new texture resolution and compute sphere radii and errors. The heights values need to be available.
	void set_texture_resolution(size_t _N)
	{
		N = _N;
		update_member(&N);
		update_tree_depth();

		errors.clear();
		radii.clear();
		errors.resize((N + 1)*(N + 1), 0.0f);
		radii.resize((N + 1)*(N + 1), 0.0f);
		processed.clear();
		processed.resize((N + 1)*(N + 1), false);

		triangle_node n = get_root_triangle(0);
		root_error = compute_error(n);
		std::cout << "root error = " << root_error << std::endl;
		std::fill(processed.begin(), processed.end(), false);
		root_radius = compute_radius(n);
		std::cout << "root radius = " << root_radius << std::endl;
		
	}
	//! function is used recursively starting from the root to compute all diamond error values 
	/*! It returns the error value of the passed node. */
	float compute_error(const triangle_node& n)
	{
		/************************************************************************************
		 tasks 3.2a: Implement a nested height error calculation on the hierarchy of right triangles
					   based on the isotropic error e_t. Do this in a recursive way over all 
					   triangles belonging to the current diamond. Store the result also in the errors-vector.
					   You can use the processed-vector to store if a certain node was already
					   calculated to avoid multiple computations of the same node.
					   Use the has_diamond_neighbor()-method to get the neighboring triangle node
		 			   of the current diamond. */

		 /*<your_code_here>*/

		//if (is_leaf(n)) return 0;
		if (level(n) >= adapted_tree_depth) return 0;
		if (processed[index(n)]) return errors[index(n)];

		float max = compute_error(left_child(n));
		max = std::max(max, compute_error(right_child(n)));

		triangle_node neighbor;
		if (has_diamond_neighbor(n, neighbor))
		{
			max = std::max(max, compute_error(left_child(neighbor)));
			max = std::max(max, compute_error(right_child(neighbor)));
		}

		short temp_length = n.omega % 2 == 0 ? n.base_length : n.base_length / 2;
		triangle_node temp_node;
		vec2 vec = nav_vertex(temp_length, (n.omega - 1) % 8);

		temp_node.x = n.x + vec[0];
		temp_node.y = n.y + vec[1];
		short h0 = heights[index(temp_node)];

		temp_node.x = n.x - vec[0];
		temp_node.y = n.y - vec[1];
		short h1 = heights[index(temp_node)];

		short h_interpolated = (h0 + h1) / 2;

		float error = std::abs(h_interpolated - heights[index(n)]);
		float result = error + max;

		errors[index(n)] = result;
		processed[index(n)] = 1;

		return result;

		/************************************************************************************/
	}
	//! function is used recursively starting from the root to compute all sphere radii
	/*! It returns the radius of the passed node. */
	float compute_radius(const triangle_node& n)
	{
		/************************************************************************************
		 tasks 3.3a: Return the radius r_t of the threshold sphere. This is also done in a
					   recursive manner as for the error computation. Store the result also in
					   the radii-vector. You can use the processed-vector to store if a certain
					   node was already calculated to avoid multiple computations of the same node.
					   Use the has_diamond_neighbor()-method to get the neighboring triangle node
		 			   of the current diamond.*/

		 /*<your_code_here>*/
		 return 0;

		/************************************************************************************/
	}
	/// based on adaptation mode check whether triangle is accurate
	bool is_accurate(const triangle_node& n) const
	{
		// if no adaption mode is defined 
		if (adaptation_mode == AM_NONE)
			return true;

		/************************************************************************************
		 tasks 3.1c: Add a check for AM_TREE_DEPTH mode: The tree level of the current node
					   has to be greater or equal than the chosen adapted_tree_depth. */

		if (level(n) < adapted_tree_depth) return false;
		 /*<your_code_here>*/

		/************************************************************************************/

		// a leaf node is accurate since there is no further refinement available
		if (is_leaf(n) || level(n) >= max_tree_depth)
			return true;

		/************************************************************************************
		 tasks 3.3c: Add criteria for the AM_ISOTROPIC_ERROR mode. Use the thresholding screen
					   space error algorithm for this which relies on the previously computed
					   errors e_t and radii r_t. The pixel error threshold is given by the UI
					   (pixel_threshold). The triangle is accurate if its pixel error is less
					   than the pixel error threshold. */

		 /*<your_code_here>*/

		/************************************************************************************/

		return true;
	}

	void subdivide(short count, vec2 base, vec2 edge0, vec2 edge1, short omega)
	{
		if (count > subdivide_count) return;
		if (count == subdivide_count)
		{
			positions.push_back(vec4(base[0], base[1], edge0[0], edge0[1]));
			normals.push_back(vec3(edge1[0], edge1[1], count));
			colors.push_back(get_orientation_color(omega));
			return;
		}

		vec2 edge0_small = vec2(edge0[0] / 2, edge0[1] / 2);
		vec2 edge1_small = vec2(edge1[0] / 2, edge1[1] / 2);

		vec2 base0 = base + edge0_small;
		vec2 base1 = base + edge1_small;
		vec2 base2 = base + edge0_small + edge1_small;

		subdivide(2 * count, base, edge0_small, edge1_small, omega);
		subdivide(2 * count, base0, edge0_small, edge1_small, omega);
		subdivide(2 * count, base1, edge0_small, edge1_small, omega);
		subdivide(2 * count, base2, -edge0_small, -edge1_small, omega);
	}

	/// extract tesselation depending on chosen adaption mode assuming diamant monotony
	void tesselate_adaptive()
	{
		/************************************************************************************
		 tasks 3.1b: Traverse through the triangle tree by using the state-less top-down refinement
					   with diamond monotonicity. Add accurate triangles to the positions, normals and
					   colors vectors that are transfered to the gpu.
					   The following informations need to be stored in these vectors:
							positions -> vec4(base.x, base.y, edge0.x, edge0.y)
   							normals   -> vec3(edge1.x, edge1.y, subdivide_count)
 							colors    -> color of triangle orientation (see get_triangle_node_color())*/

		 /*<your_code_here>*/
		std::queue<triangle_node> queue;
		queue.push(get_root_triangle(0));
		queue.push(get_root_triangle(1));

		while (!queue.empty())
		{
			triangle_node node = queue.front();
			queue.pop();

			if (node.x < 0) node.x = 0;
			if (node.x > N) node.x = N;
			if (node.y < 0) node.y = 0;
			if (node.y > N) node.y = N;

			if (is_accurate(node))
			{
				vec2 base;
				vec2 edge0;
				vec2 edge1;

				short sgn_x, sgn_y;

				if (node.omega % 2 == 1)
				{
					sgn_x = node.omega == 3 ? 1 : node.omega == 7 ? -1 : 0;
					sgn_y = node.omega == 1 ? -1 : node.omega == 5 ? 1 : 0;
				}
				else
				{
					sgn_x = node.omega == 0 || node.omega == 6 ? -1 : 1;
					sgn_y = node.omega == 0 || node.omega == 2 ? -1 : 1;
				}

				short distance_to_side = node.base_length / 2;

				base[0] = node.x + sgn_x * distance_to_side;
				base[1] = node.y + sgn_y * distance_to_side;
				edge0 = nav_vertex(node.base_length, node.omega);
				edge1 = nav_vertex(node.base_length, (node.omega + 2) % 8);

				subdivide(1, base, edge0, edge1, node.omega);

				/************************************************************************************
			   tasks 3.2b: If show_error_spheres is true fill the spheres vector with the world point
			   of accurate triangle nodes. Set the error multiplied with the z-extent as
			   a radius and fill the spheres_colors vector with the color of the triangle
			   node.
			   The following informations need to be stored in the spheres-vector:
			   vec4(n.x, n.y, n.z, radius)*/

			   /*<your_code_here>*/
				//compute_error(node);
				vec3 point = world_point(node);
				float radius = errors[index(node)] * extent(2);
				spheres.push_back(vec4(point, radius));
				sphere_colors.push_back(get_orientation_color(node.omega));

			  /************************************************************************************/
			}
			else
			{
				queue.push(left_child(node));
				queue.push(right_child(node));
			}
		}

		/************************************************************************************/

		/************************************************************************************
		 tasks 3.3b: If show_radius_spheres is true fill the spheres vector with the world point of
					   the accurate triangle nodes and set the radius as the sphere radius.
					   Fill the spheres_colors vector with the color of the triangle node.
					   The following informations need to be stored in the spheres-vector:
					   vec4(n.x, n.y, n.z, radius)*/

		 /*<your_code_here>*/

		/************************************************************************************/
	}
	/// resulting number of triangles
	int nr_triangles;
	//@}

	/**@name file io*/
	//@{
private:
	/// remember whether to read the dem texture in the init_frame method
	bool read_dem_file;
	/// remember whether to read the color texture in the init_frame method
	bool read_color_file;
protected:
	/// file name of dem texture
	std::string dem_file_name;
	/// file name of color texture
	std::string color_file_name;
	/// read a texture from a given file into given texture object and extract height array in case of dem image
	bool read_texture(cgv::render::context& ctx, const std::string& file_name, cgv::render::texture& tex, bool extract_heights = false)
	{
		cgv::data::data_format fmt;
		cgv::data::data_view dv;
		cgv::media::image::image_reader ir(fmt);
		if (!ir.open(file_name)) {
			std::cerr << "could not open file " << file_name << std::endl;
			return false;
		}
		if (!ir.read_image(dv)) {
			std::cerr << "could not read file " << file_name << std::endl;
			return false;
		}
		if (tex.is_created())
			tex.destruct(ctx);
		tex.create(ctx, dv);
		if (extract_heights) {
			size_t N = tex.get_width();
			size_t idx = 0;
			heights.resize((N + 1)*(N + 1));
			if (fmt.get_component_type() == cgv::type::info::TI_UINT16) {
				max_dem_value = 65535;
				unsigned short* ptr = dv.get_ptr<unsigned short>();
				unsigned incr = fmt.get_entry_size() / 2;
				for (size_t j = 0; j < N; ++j) {
					for (size_t i = 0; i < N; ++i) {
						heights[idx] = *ptr;
						++idx;
						ptr += incr;
					}
					heights[idx] = heights[idx - 1];
					++idx;
				}
			}
			else {
				max_dem_value = 255;
				unsigned char* ptr = dv.get_ptr<unsigned char>();
				unsigned incr = fmt.get_entry_size();
				for (size_t j = 0; j < N; ++j) {
					for (size_t i = 0; i < N; ++i) {
						heights[idx] = *ptr;
						++idx;
						ptr += incr;
					}
					heights[idx] = heights[idx - 1];
					++idx;
				}
			}
			for (size_t i = 0; i < N; ++i) {
				heights[idx] = heights[idx - N - 1];
				++idx;
			}
			heights[idx] = heights[idx - 1];
		}
		std::cout << "read " << file_name << " (" << fmt.get_width() << "x" << fmt.get_height() << ")" << std::endl;
		return true;
	}
	//@}

	/**@name rendering*/
	//@{
	/// pointer to view
	cgv::render::view* view_ptr;
	/// shader program
	cgv::render::shader_program terrain_prog;
	/// whether to construct triangular batches with instanced rendering (normally true)
	bool triangular;
	/// dem texture
	cgv::render::texture dem_tex;
	/// dem filter parameters
	cgv::render::TextureFilter dem_minification, dem_magnification;
	/// color texture
	cgv::render::texture color_tex;
	/// color filter parameters
	cgv::render::TextureFilter color_minification, color_magnification;
	/// position vertex attribute container
	std::vector<vec4> positions;
	/// normal vertex attribute container
	std::vector<vec3> normals;
	/// color vertex attribute container
	std::vector<rgb>  colors;
	/// set the global attribute array bindings to the vertex attribute containers
	void set_vertex_attributes(cgv::render::context& ctx)
	{
		if (positions.size() > 0) {
			int pos_idx = terrain_prog.get_attribute_location(ctx, "position");
			cgv::render::attribute_array_binding::set_global_attribute_array(ctx, pos_idx, positions);
			cgv::render::attribute_array_binding::enable_global_array(ctx, pos_idx);
		}
		if (normals.size() > 0) {
			int nml_idx = terrain_prog.get_attribute_location(ctx, "normal");
			cgv::render::attribute_array_binding::set_global_attribute_array(ctx, nml_idx, normals);
			cgv::render::attribute_array_binding::enable_global_array(ctx, nml_idx);
		}
		if (colors.size() > 0) {
			int clr_idx = terrain_prog.get_attribute_location(ctx, "color");
			cgv::render::attribute_array_binding::set_global_attribute_array(ctx, clr_idx, colors);
			cgv::render::attribute_array_binding::enable_global_array(ctx, clr_idx);
		}
	}
	/// vertex attribute container for spheres
	std::vector<vec4> spheres;
	/// vertex attribute container for sphere colors
	std::vector<rgb>  sphere_colors;
	
	//@}

	/**@name visualization*/
	//@{
	/// blend in error based color mapping
	float error_lambda;
	/// blend in radius based color mapping
	float radius_lambda;
	/// compute triangle color based on the lambdas for color, error and radius 
	static rgb get_orientation_color(short omega) 
	{
		static rgb colors[8] = {
			{ 1.0f, 0.0f, 0.0f },
			{ 0.5f, 0.5f, 0.0f },
			{ 1.0f, 1.0f, 0.0f },
			{ 0.5f, 1.0f, 0.0f },
			{ 0.0f, 1.0f, 0.0f },
			{ 0.0f, 0.5f, 0.5f },
			{ 0.0f, 0.0f, 1.0f },
			{ 0.5f, 0.0f, 0.5f }
		};
		return colors[omega];
	}
	/// compute triangle color based on the lambdas for color, error and radius 
	rgb get_triangle_node_color(const triangle_node& n)
	{
		unsigned idx = n.y*(N + 1) + n.x;
		float error = errors[idx];
		float v = error / root_error;
		float radius = radii[idx];
		float w = radius / root_radius;
		rgb error_color(v, 0.5f, 1.0f - v);
		rgb radius_color(0.5f, w, 0.5f);
		float lambda = std::max(error_lambda, radius_lambda);
		if (lambda < 0.001f)
			return get_orientation_color(n.omega);
		rgb color = error_lambda / (error_lambda + radius_lambda)*error_color +
			radius_lambda / (error_lambda + radius_lambda)*radius_color;
		return lambda * color + (1 - lambda)*get_orientation_color(n.omega);
	}
	/// lambda to blend between color from color texture and from error/radius/orientation visualization
	float color_lambda;
	/// whether to additionally show wireframe
	bool wireframe;
	/// threshold used to control width of wireframe
	float wire_threshold;
	/// color of wireframe
	rgb wire_color;
	/// whether to show radius spheres
	bool show_radius_spheres;
	/// whether to show error spheres
	bool show_error_spheres;
	/// style for rendering spheres
	cgv::render::sphere_render_style srs;
	//@}

public:
	terrain() :
		node("terrain")
	{
		max_tree_depth = 48;
		adapted_tree_depth = 5;
		adaptation_mode = AM_TREE_DEPTH;
		root_error = 0;
		root_radius = 0;
		pixel_threshold = 5;
		max_dem_value = 65535;

		view_ptr = 0;
		view_factor = 1;
		eye_world = vec3(0.0f);

		color_lambda = 0.5f;
		error_lambda = 0.5f;
		radius_lambda = 0.5;

		extent = vec3(20, 10, 0.04f);
		triangular = true;
		subdivide_count = 1;
		wireframe = false;
		wire_threshold = 0.02f;
		wire_color = rgb(0.2f,0.2f,0.2f);
		show_radius_spheres = false;
		show_error_spheres = false;

		read_dem_file = false;
		read_color_file = false;
		dem_tex.set_min_filter  (dem_minification  = cgv::render::TF_LINEAR_MIPMAP_LINEAR);
		dem_tex.set_mag_filter(dem_magnification   = cgv::render::TF_LINEAR);
		dem_tex.set_wrap_r(cgv::render::TW_CLAMP_TO_EDGE);
		dem_tex.set_wrap_s(cgv::render::TW_CLAMP_TO_EDGE);
		color_tex.set_min_filter(color_minification  = cgv::render::TF_LINEAR_MIPMAP_LINEAR);
		color_tex.set_mag_filter(color_magnification = cgv::render::TF_LINEAR);
		color_tex.set_wrap_r(cgv::render::TW_CLAMP_TO_EDGE);
		color_tex.set_wrap_s(cgv::render::TW_CLAMP_TO_EDGE);
	}
	/// return type name
	std::string get_type_name() const
	{
		return "terrain";
	}
	/// makes data member of this class available as named properties (which can be set in a config file)
	bool self_reflect(cgv::reflect::reflection_handler& rh)
	{
		return
			rh.reflect_member("dem_file_name", dem_file_name) &&
			rh.reflect_member("show_radius_spheres", show_radius_spheres) &&
			rh.reflect_member("show_radius_spheres", show_error_spheres) &&
			rh.reflect_member("color_file_name", color_file_name) &&
			rh.reflect_member("subdivide_count", subdivide_count);
	}
	/// callback for all changed UI elements
	void on_set(void* member_ptr)
	{
		if (member_ptr == &dem_file_name) {
			read_dem_file = true;
		}
		if (member_ptr >= &extent && member_ptr < &extent+1) {
			std::fill(processed.begin(), processed.end(), false);
			root_radius = compute_radius(get_root_triangle(0));
			std::cout << "root radius = " << root_radius << std::endl;
			radii_cache = radii;
		}

		if (member_ptr == &color_file_name) {
			read_color_file = true;
		}
		if (member_ptr == &subdivide_count) {
			update_tree_depth();
		}
		if (member_ptr == &dem_minification)
			dem_tex.set_min_filter(dem_minification);
		if (member_ptr == &dem_magnification)
			dem_tex.set_mag_filter(dem_magnification);
		if (member_ptr == &color_minification)
			color_tex.set_min_filter(color_minification);
		if (member_ptr == &color_magnification)
			color_tex.set_mag_filter(color_magnification);

		update_member(member_ptr);
		post_redraw();
	}
	/// initialize everything that needs the context
	bool init(cgv::render::context& ctx)
	{
		view_ptr = find_view_as_node();
		// set background color to white
		ctx.set_bg_clr_idx(3);
		// build shader program for spheres
		if (!terrain_prog.build_program(ctx, "terrain.glpr")) {
			std::cerr << "could not build terrain shader program" << std::endl;
			exit(0);
		}
		cgv::render::ref_sphere_renderer(ctx, 1);
		return true;
	}
	/// destruct all graphics objects
	void clear(cgv::render::context& ctx)
	{
		dem_tex.destruct(ctx);
		color_tex.destruct(ctx);
		terrain_prog.destruct(ctx);
		cgv::render::ref_sphere_renderer(ctx, -1);
	}
	/// read textures here to have access to the context
	void init_frame(cgv::render::context& ctx)
	{
		// ensure that textures are read
		if (read_dem_file) {
			read_texture(ctx, dem_file_name, dem_tex, true);
			set_texture_resolution(dem_tex.get_width());
			read_dem_file = false;
		}
		if (read_color_file) {
			read_texture(ctx, color_file_name, color_tex);
			read_color_file = false;
		}
	}
	/// this method is called to draw a frame
	void draw(cgv::render::context& ctx)
	{
		if (!terrain_prog.is_linked())
			return;

		eye_world = view_ptr->get_eye();
		view_factor = float((2.0f*view_ptr->get_tan_of_half_of_fovy(false))/ctx.get_height());

		// change nothing if adaptive mode is set to AM_NONE
		if (adaptation_mode != AM_NONE || positions.empty()) {
			positions.clear();
			normals.clear();
			colors.clear();
			spheres.clear();
			sphere_colors.clear();
			tesselate_adaptive();
			nr_triangles = positions.size()*subdivide_count*subdivide_count;
			update_member(&nr_triangles);
		}
		set_vertex_attributes(ctx);

		// render terrain
		terrain_prog.enable(ctx);
		if (dem_tex.is_created())
			dem_tex.enable(ctx, 0);
		if (color_tex.is_created())
			color_tex.enable(ctx, 1);
			vec2 scale = vec2(1.0f / float(dem_tex.get_width()), 1.0f / float(dem_tex.get_height()));
			terrain_prog.set_uniform(ctx, "dem_tex", 0);
			terrain_prog.set_uniform(ctx, "triangular", triangular);
			terrain_prog.set_uniform(ctx, "wireframe", wireframe);
			terrain_prog.set_uniform(ctx, "wire_threshold", wire_threshold);
			terrain_prog.set_uniform(ctx, "wire_color", wire_color);
			terrain_prog.set_uniform(ctx, "map_color_to_material", 3);
			terrain_prog.set_uniform(ctx, "color_tex", 1);
			terrain_prog.set_uniform(ctx, "color_lambda", color_lambda);
			terrain_prog.set_uniform(ctx, "scale", scale);
			terrain_prog.set_uniform(ctx, "extent", extent);
			terrain_prog.set_uniform(ctx, "N", (int)N);
			glDrawArraysInstanced(GL_POINTS, 0, positions.size(), subdivide_count*subdivide_count);
		if (color_tex.is_created())
			color_tex.disable(ctx);
		if (dem_tex.is_created())
			dem_tex.disable(ctx);
		terrain_prog.disable(ctx);

		// render spheres
		if (spheres.size() > 0) {
			cgv::render::sphere_renderer& sr = cgv::render::ref_sphere_renderer(ctx);
			sr.set_render_style(srs);
			sr.set_sphere_array(ctx, spheres);
			sr.set_color_array(ctx, sphere_colors);
			if (sr.validate_and_enable(ctx)) {
				glDrawArrays(GL_POINTS, 0, (GLsizei)spheres.size());
				sr.disable(ctx);
			}
		}
	}
	/// not yet implemented
	bool handle(cgv::gui::event& e)
	{
		if (e.get_kind() != cgv::gui::EID_KEY)
			return false;
		cgv::gui::key_event& ke = static_cast<cgv::gui::key_event&>(e);
		if (ke.get_action() == cgv::gui::KA_RELEASE)
			return false;
		switch (ke.get_key()) {
		case 'A':
			if (ke.get_modifiers() == 0) {
				if (adaptation_mode == AM_LAST)
					adaptation_mode = AM_FIRST;
				else
					++(int&)adaptation_mode;
				on_set(&adaptation_mode);
				return true;
			}
			else if (ke.get_modifiers() == cgv::gui::EM_SHIFT) {
				if (adaptation_mode == AM_FIRST)
					adaptation_mode = AM_LAST;
				else
					--(int&)adaptation_mode;
				on_set(&adaptation_mode);
				return true;
			}
			break;
		case 'C':
			if (ke.get_modifiers() == 0) {
				srs.culling_mode = srs.culling_mode == cgv::render::CM_OFF ?
					cgv::render::CM_FRONTFACE : cgv::render::CM_OFF;
				on_set(&srs.culling_mode);
				return true;
			}
			break;
		case 'D':
			if (ke.get_modifiers() == 0) {
				if (adapted_tree_depth < tree_depth) {
					++adapted_tree_depth;
					on_set(&adapted_tree_depth);
				}
				return true;
			}
			else if (ke.get_modifiers() == cgv::gui::EM_SHIFT) {
				if (adapted_tree_depth > 0) {
					--adapted_tree_depth;
					on_set(&adapted_tree_depth);
				}
				return true;
			}
			break;
		case 'S':
			if (ke.get_modifiers() == 0) {
				if (subdivide_count < 256) {
					subdivide_count *= 2;
					on_set(&subdivide_count);
				}
				return true;
			}
			else if (ke.get_modifiers() == cgv::gui::EM_SHIFT) {
				if (subdivide_count > 1) {
					subdivide_count /= 2;
					on_set(&subdivide_count);
				}
				return true;
			}
			break;
		case 'R':
			if (ke.get_modifiers() == 0) {
				show_radius_spheres = !show_radius_spheres;
				on_set(&show_radius_spheres);
				return true;
			}
			break;
		case 'W':
			if (ke.get_modifiers() == 0) {
				wireframe = !wireframe;
				on_set(&wireframe);
				return true;
			}
			break;
		case 'E':
			if (ke.get_modifiers() == 0) {
				show_error_spheres = !show_error_spheres;
				on_set(&show_error_spheres);
				return true;
			}
			break;
		}
		return false;
	}
	/// not yet implemented
	void stream_help(std::ostream& os)
	{
		os << "terrain: ..." << std::endl;
	}
	/// not yet implemented
	void stream_stats(std::ostream& os)
	{

	}
	/// 
	void create_gui()
	{
		add_decorator("terrain", "heading");

		bool parameter = true;
		if (begin_tree_node("io", parameter, false)) {
			align("\a");
			add_gui("dem_file", dem_file_name, "file_name", "title='open dem image';filter='images (bmp,jpg,png,tif):*.bmp,*.jpg,*.png,*.tif|all files:*.*'");
			add_gui("color_file", color_file_name, "file_name", "title='open color image';filter='images (bmp,jpg,png,tif):*.bmp,*.jpg,*.png,*.tif|all files:*.*'");
			align("\b");
			end_tree_node(parameter);
		}
		if (begin_tree_node("adaptive hierarchy", max_tree_depth, true)) {
			align("\a");
			add_member_control(this, "extent_x", extent[0], "value_slider", "ticks=true;min=1;max=100;log=true");
			add_member_control(this, "extent_y", extent[1], "value_slider", "ticks=true;min=1;max=100;log=true");
			add_member_control(this, "extent_z", extent[2], "value_slider", "ticks=true;min=0.01;max=10;log=true");
			add_view("N", N);
			add_view("nr_triangles", nr_triangles);
			add_member_control(this, "subdivide_count", (cgv::type::DummyEnum&)subdivide_count, "dropdown", "enums='1=1,2=2,4=4,8=8,16=16,32=32,64=64,128=128");
			add_member_control(this, "adaptation_mode", adaptation_mode, "dropdown", "enums='none,tree_depth,isotropic,anisotropic'");
			add_member_control(this, "adapted_tree_depth", adapted_tree_depth, "value_slider", "min=0;max=12");
			add_member_control(this, "pixel_threshold", pixel_threshold, "value_slider", "min=0;max=20;log=true;ticks=true");
			add_member_control(this, "max_tree_depth", max_tree_depth, "value_slider", "min=1;max=12");
			configure_gui();
			align("\b");
			end_tree_node(max_tree_depth);
		}
		if (begin_tree_node("rendering", extent, false)) {
			align("\a");
			add_member_control(this, "triangular", triangular, "toggle");
			add_member_control(this, "dem_minification", dem_minification, "dropdown",
				"enums='nearest,linear,nearest_mm_nearest,linear_mm_nearest,linear_mm_linear,anisotrop'");
			add_member_control(this, "dem_magnification", dem_magnification, "dropdown", "enums='nearest,linear'");
			add_member_control(this, "color_minification", color_minification, "dropdown",
				"enums='nearest,linear,nearest_mm_nearest,linear_mm_nearest,linear_mm_linear,anisotrop'");
			add_member_control(this, "color_magnification", color_magnification, "dropdown", "enums='nearest,linear'");
			align("\b");
			end_tree_node(extent);
		}
		if (begin_tree_node("visualization", color_lambda, true)) {
			align("\a");
			add_member_control(this, "color_lambda", color_lambda, "value_slider", "ticks=true;min=0;max=1");
			add_member_control(this, "error_lambda", error_lambda, "value_slider", "ticks=true;min=0;max=1");
			add_member_control(this, "radius_lambda", radius_lambda, "value_slider", "ticks=true;min=0;max=1");
			add_member_control(this, "wireframe", wireframe, "toggle");
			add_member_control(this, "wire_threshold", wire_threshold, "value_slider", "ticks=true;min=0;max=1;log=true");
			add_member_control(this, "wire_color", wire_color);
			add_member_control(this, "show_error_spheres", show_error_spheres, "toggle");
			add_member_control(this, "show_radius_spheres", show_radius_spheres, "toggle");
			if (begin_tree_node("sphere style", srs, false)) {
				align("\a");
				add_gui("srs", srs);
				align("\b");
				end_tree_node(color_lambda);
			}
			align("\b");
			end_tree_node(color_lambda);
		}
	}
};

cgv::base::object_registration<terrain> terrain_reg("terrain");

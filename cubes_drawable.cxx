// Framework core
#include <cgv/base/register.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/trigger.h>
#include <cgv/render/drawable.h>
#include <cgv/render/render_types.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/texture.h>
#include <cgv/render/frame_buffer.h>
#include <cgv/render/vertex_buffer.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv/media/font/font.h>
#include <cgv/math/ftransform.h>

// Framework standard plugins
#include <cgv_gl/gl/gl.h>

#include <cubes_fractal.h>

const int VERTEX_COUNT = 14,
VERTEX_ARRAY_STRIDE = 3;

class cgv_cubes
	: public cgv::base::base,      // This class supports reflection
	public cgv::gui::provider,   // Instances of this class provde a GUI
	public cgv::render::drawable // Instances of this class can be rendered
{

protected:



	unsigned int recursion_depth;
	enum VaMode { NO_VERTEX_ARRAY, INTERLEAVED, NON_INTERLEAVED } va_mode;
	float cube_r, cube_g, cube_b;

	cgv::media::color<float> cube_color;
	cubes_fractal fractal;

	struct vertex {
		cgv::render::render_types::vec3 pos;
		cgv::render::render_types::vec2 tcoord;
	};
	std::vector<vertex> vertices;
	cgv::render::vertex_buffer vb;
	cgv::render::attribute_array_binding vertex_array;



public:

	cgv_cubes() : cube_r(0.0f), cube_g(0.0f), cube_b(0.0f), cube_color(cube_r, cube_g, cube_b)
	{}

	// Should be overwritten to sensibly implement the cgv::base::named interface
	std::string get_type_name(void) const
	{
		return "cgv_cubes";
	}

	// Part of the cgv::base::base interface, can be implemented to make data members of
	// this class available as named properties, e.g. for use with config files
	bool self_reflect(cgv::reflect::reflection_handler& rh)
	{
		unsigned *va_mode_uint = (unsigned*)&va_mode;

		// Reflect the properties
		return
			rh.reflect_member("cube_r", cube_r) &&
			rh.reflect_member("cube_g", cube_g) &&
			rh.reflect_member("cube_b", cube_b) &&
			rh.reflect_member("recursion_depth", recursion_depth) &&
			rh.reflect_member("va_mode", *va_mode_uint);

	}

	// Part of the cgv::base::base interface, should be implemented to respond to write
	// access to reflected data members of this class, e.g. from config file processing
	// or gui interaction.
	void on_set(void* member_ptr)
	{
		if (member_ptr == &cube_r || member_ptr == &cube_g ||
			member_ptr == &cube_b)
		{
			cube_color.R() = cube_r;
			cube_color.G() = cube_g;
			cube_color.B() = cube_b;
			update_member(&cube_color);
		}
		// ...and vice versa (can only happen via GUI interaction)
		if (member_ptr == &cube_color)
		{
			cube_r = cube_color.R();
			cube_g = cube_color.G();
			cube_b = cube_color.B();
		}

		if (member_ptr == &va_mode) {
			if (va_mode == NO_VERTEX_ARRAY) fractal.use_vertex_array(nullptr, 0, GL_QUADS);
			else fractal.use_vertex_array(&vertex_array, vertices.size(), GL_TRIANGLE_STRIP);
		}


		// Make sure the GUI reflects the new state, in case the write access did not
		// originate from GUI interaction
		update_member(member_ptr);

		// Also trigger a redraw in case the drawable node is active
		if (this->is_visible())
			post_redraw();
	}

	// We use this for validating GUI input
	bool gui_check_value(cgv::gui::control<int>& ctrl)
	{

		// Check passed
		return true;
	}

	// We use this for acting upon validated GUI input
	void gui_value_changed(cgv::gui::control<int>& ctrl)
	{

		// Redraw the scene
		post_redraw();
	}

	// Required interface for cgv::gui::provider
	void create_gui(void)
	{
		// Simple controls. Notifies us of GUI input via the on_set() method.
		// - section header
		add_decorator("Cube properties", "heading", "level=1");

		add_member_control(this, "recursion depth", recursion_depth);
		add_member_control(this, "cube color", cube_color);
		add_member_control(
			this, "vertex array", va_mode, "dropdown",
			"enums='No Vertex Array,Interleaved,Non-Interleaved'"
		);
	}

	// Part of the cgv::render::drawable interface, can be overwritten if there is some
	// intialization work to be done that needs a set-up and ready graphics context,
	// which usually you don't have at object construction time. Should return true if
	// the initialization was successful, false otherwise.
	bool init(cgv::render::context& ctx)
	{
		// Keep track of success - do it this way (instead of e.g. returning false
		// immediatly) to perform every init step even if some go wrong.
		bool success = true;

		init_unit_cube_geometry();

		cgv::render::shader_program& default_shader
			= ctx.ref_default_shader_program(true);

		cgv::render::type_descriptor
			vec2type =
			cgv::render::element_descriptor_traits<cgv::render::render_types::vec2>
			::get_type_descriptor(vertices[0].tcoord),
			vec3type =
			cgv::render::element_descriptor_traits<cgv::render::render_types::vec3>
			::get_type_descriptor(vertices[0].pos);
		// - create buffer objects
		success = vb.create(ctx, &(vertices[0]), vertices.size()) && success;
		success = vertex_array.create(ctx) && success;
		success = vertex_array.set_attribute_array(
			ctx, default_shader.get_position_index(), vec3type, vb,
			0, // position is at start of the struct <-> offset = 0
			vertices.size(), // number of position elements in the array
			sizeof(vertex) // stride from one element to next
		) && success;
		success = vertex_array.set_attribute_array(
			ctx, default_shader.get_texcoord_index(), vec2type, vb,
			sizeof(vertex::pos), // texture coords come after position in our struct
			vertices.size(), // number of texcoord elements in the array
			sizeof(vertex) // stride from one element to next
		) && success;


		// All initialization has been attempted
		return success;
	}

	void draw(cgv::render::context& ctx)
	{
		cgv::render::shader_program& default_shader = ctx.ref_surface_shader_program(false);
		default_shader.enable(ctx);
		ctx.set_color(rgb(1.0f));

		/* if (use_vertex_array) {

			unsigned int num_vertices;
			switch (recursion_depth) {
			case 0:
				num_vertices = 8 * 1;
				break;
			case 1:
				num_vertices = 8 * (1 + 4);
				break;

			default:
				num_vertices = 8 * (1 + 4 + 3 * (recursion_depth - 2));
				break;
			}

			bool success = true;
		}*/

		fractal.draw_recursive(ctx, cube_color, recursion_depth, 0);

		default_shader.disable(ctx);
	}

	void init_unit_cube_geometry(void)
	{


		// http://www.cs.umd.edu/gvil/papers/av_ts.pdf
		// https://stackoverflow.com/questions/28375338/cube-using-single-gl-triangle-strip
		// Efficient GL_TRIANGLE_STRIP method
		static float vertices_data_array[VERTEX_COUNT * VERTEX_ARRAY_STRIDE] = {

			// FRONT
			-1.f, 1.f, 1.f,     // Front-top-left
			1.f, 1.f, 1.f,      // Front-top-right
			-1.f, -1.f, 1.f,    // Front-bottom-left
			1.f, -1.f, 1.f,     // Front-bottom-right
			1.f, -1.f, -1.f,    // Back-bottom-right
			1.f, 1.f, 1.f,      // Front-top-right
			1.f, 1.f, -1.f,     // Back-top-right
			-1.f, 1.f, 1.f,     // Front-top-left
			-1.f, 1.f, -1.f,    // Back-top-left
			-1.f, -1.f, 1.f,    // Front-bottom-left
			-1.f, -1.f, -1.f,   // Back-bottom-left
			1.f, -1.f, -1.f,    // Back-bottom-right
			-1.f, 1.f, -1.f,    // Back-top-left
			1.f, 1.f, -1.f      // Back-top-right

		};

		vertices.resize(VERTEX_COUNT);
		for (int i = 0; i < VERTEX_COUNT; i += 1) {
			int current_data_array_index_base = i * VERTEX_ARRAY_STRIDE;
			vertices[i].pos.set(
				vertices_data_array[current_data_array_index_base],
				vertices_data_array[current_data_array_index_base + 1],
				vertices_data_array[current_data_array_index_base + 2]
			);
			vertices[i].tcoord.set(
				i / VERTEX_COUNT,
				i / VERTEX_COUNT
			);
		}
	}


};

// Create an instance of the cubes class at plugin load and register it with the framework
cgv::base::object_registration<cgv_cubes> cgv_cubes_registration("");
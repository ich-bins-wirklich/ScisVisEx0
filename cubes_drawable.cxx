// This source code is property of the Computer Graphics and Visualization chair of the
// TU Dresden. Do not distribute! 
// Copyright (C) CGV TU Dresden - All Rights Reserved
//
// The main file of the plugin. It defines a class that demonstrates how to register with
// the scene graph, drawing primitives, creating a GUI, using a config file and various
// other parts of the framework.

// Framework core
#include <cgv/base/register.h>
#include <cgv/gui/provider.h>
#include <cgv/gui/trigger.h>
#include <cgv/render/drawable.h>
#include <cgv/render/render_types.h>
#include <cgv/render/shader_program.h>
#include <cgv/render/vertex_buffer.h>
#include <cgv/render/attribute_array_binding.h>
#include <cgv/math/ftransform.h>

// Framework standard plugins
#include <cgv_gl/gl/gl.h>

// Local includes
#include "cubes_fractal.h"



/// ************************************************************************************/
/// Task 1.2a: Create a drawable that provides a (for now, empty) GUI and supports
///            reflection, so that its properties can be set via config file.
///
/// Task 1.2b: Utilize the cubes_fractal class to render a fractal of hierarchically
///            transformed cubes. Expose its recursion depth and color properties to GUI
///            manipulation and reflection. Set reasonable values via the config
///            file.
///
/// Task 1.2c: Implement an option (configurable via GUI and config file) to use a vertex
///            array object for rendering the cubes. The vertex array functionality 
///            should support (again, configurable via GUI and config file) both
///            interleaved (as in cgv_demo.cpp) and non-interleaved attributes.

class cgv_cubes
	: public cgv::base::base,      // This class supports reflection
	public cgv::gui::provider,   // Instances of this class provde a GUI
	public cgv::render::drawable // Instances of this class can be rendered
{

protected:
	unsigned int recursion_depth;
	float cube_r, cube_g, cube_b;
	enum VaMode { NO_VERTEX_ARRAY, INTERLEAVED, NON_INTERLEAVED } va_mode;

	cgv::media::color<float> cube_color;
	cubes_fractal fractal;


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


		// Make sure the GUI reflects the new state, in case the write access did not
		// originate from GUI interaction
		update_member(member_ptr);

		// Also trigger a redraw in case the drawable node is active
		if (this->is_visible())
			post_redraw();
	}

	// We use this for validating GUI input
	bool gui_check_value(cgv::gui::control<int> &ctrl)
	{

		// Check passed
		return true;
	}

	// We use this for acting upon validated GUI input
	void gui_value_changed(cgv::gui::control<int> &ctrl)
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

		// All initialization has been attempted
		return success;
	}

	void draw(cgv::render::context& ctx)
	{
		cgv::render::shader_program &default_shader = ctx.ref_surface_shader_program(false);
		default_shader.enable(ctx);

		fractal.draw_recursive(ctx, cube_color, recursion_depth, 0);

		default_shader.disable(ctx);
	}


};

/// [END] Tasks 1.2a, 1.2b and 1.2c
/// ************************************************************************************/



/// ************************************************************************************/
/// Task 1.2a: register an instance of your drawable.
// Create an instance of the cubes class at plugin load and register it with the framework
cgv::base::object_registration<cgv_cubes> cgv_cubes_registration("");
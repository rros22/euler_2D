#include <stdio.h>
#include <math.h>
#include "godunov/godunov.h"
#include "post_processing/post_processing.h"
#include "box/box.h"

int main()
{
    const int base_x = 100;
    const int base_y = 100;

    const int nx_side = 2 * base_x;
    const int nx_bottom = 1 * base_x;
    const int ny_bottom = 6 * base_y;
    const int nx_top = nx_bottom;
    const int ny_top = 2 * base_y;

    const double Lx = 1.0;
    const double Ly = 1.0;
    const double dx = Lx / (nx_side + nx_bottom + nx_side); // Arbitrary scale
    const double dy = Ly / (ny_top + ny_bottom);

    const double cfl = 0.4;
    const double gamma = 1.4;

    // Initial conditions
    double w_left[4] = {0.125, 0.0, 0.0, 0.1};
    double w_right[4] = {0.125, 0.0, 0.0, 0.1};
    double w_bottom[4] = {0.125, 0.0, 0.0, 0.1};
    double w_top[4] = {1.0, 0.0, 0.0, 1.0};

    // Create boxes
    struct Box *box_left = box_create(nx_side, ny_bottom);
    struct Box *box_bottom = box_create(nx_bottom, ny_bottom);
    struct Box *box_right = box_create(nx_side, ny_bottom);
    struct Box *box_top = box_create(nx_top, ny_top);

    struct Box *temp_left = box_create(nx_side, ny_bottom);
    struct Box *temp_bottom = box_create(nx_bottom, ny_bottom);
    struct Box *temp_right = box_create(nx_side, ny_bottom);
    struct Box *temp_top = box_create(nx_top, ny_top);

    struct Box *boxes[] = {box_left, box_bottom, box_right, box_top};
    struct Box *temps[] = {temp_left, temp_bottom, temp_right, temp_top};
    const int num_boxes = 4;

    // Connections
    box_set_boundary_connected(box_left, EDGE_RIGHT, box_bottom, EDGE_LEFT);
    box_set_boundary_connected(box_bottom, EDGE_LEFT, box_left, EDGE_RIGHT);

    box_set_boundary_connected(box_bottom, EDGE_RIGHT, box_right, EDGE_LEFT);
    box_set_boundary_connected(box_right, EDGE_LEFT, box_bottom, EDGE_RIGHT);

    box_set_boundary_connected(box_bottom, EDGE_TOP, box_top, EDGE_BOTTOM);
    box_set_boundary_connected(box_top, EDGE_BOTTOM, box_bottom, EDGE_TOP);

    // Transmissive boundaries
    box_set_boundary_transmissive(box_left, EDGE_LEFT);
    box_set_boundary_transmissive(box_right, EDGE_RIGHT);
    box_set_boundary_transmissive(box_left, EDGE_BOTTOM);
    box_set_boundary_transmissive(box_right, EDGE_BOTTOM);
    box_set_boundary_transmissive(box_left, EDGE_TOP);
    box_set_boundary_transmissive(box_right, EDGE_TOP);
    box_set_boundary_transmissive(box_bottom, EDGE_BOTTOM);

    // Reflective boundaries
    box_set_boundary_reflective(box_top, EDGE_TOP);
    box_set_boundary_reflective(box_top, EDGE_LEFT);
    box_set_boundary_reflective(box_top, EDGE_RIGHT);

    // Copy boundaries
    for (int i = 0; i < num_boxes; ++i)
        box_copy_boundaries(temps[i], boxes[i]);

    // Set initial conditions for each box
    box_set_uniform_initial_conditions(box_left, w_left);
    box_set_uniform_initial_conditions(box_right, w_right);
    box_set_uniform_initial_conditions(box_bottom, w_bottom);
    box_set_uniform_initial_conditions(box_top, w_top);
    // Convert to conserved
    for (int i = 0; i < num_boxes; ++i)
        box_prim_to_cons(boxes[i], gamma);

    // Save initial
    save_box_to_csv(box_left, "results/box_left_initial.csv");
    save_box_to_csv(box_bottom, "results/box_bottom_initial.csv");
    save_box_to_csv(box_right, "results/box_right_initial.csv");
    save_box_to_csv(box_top, "results/box_top_initial.csv");
    save_box_connections(boxes, num_boxes, "results/box_connections.csv");

    // Main loop
    double t = 0.0;
    int step = 0;
    const int max_steps = 20000;
    const double stop_time = 10;

    while (step < max_steps && t < stop_time)
    {
        printf("Step %d: Begin\n", step);

        double vmax = 0.0;
        for (int i = 0; i < num_boxes; ++i)
        {
            double v = box_max_wave_speed(boxes[i], gamma);
            if (v > vmax)
                vmax = v;
        }

        double dt = cfl * fmin(dx, dy) / vmax;
        printf("  Time step dt = %.5e\n", dt);

        // X-sweep
        for (int i = 0; i < num_boxes; ++i)
            box_apply_boundary_x_edges(boxes[i]);
        for (int i = 0; i < num_boxes; ++i)
            godunov_x_sweep(boxes[i], temps[i], dt, dx, gamma);
        for (int i = 0; i < num_boxes; ++i)
        {
            box_cons_to_prim(temps[i], gamma);
            box_swap(&boxes[i], &temps[i]);
        }

        // Y-sweep 1
        for (int i = 0; i < num_boxes; ++i)
            box_apply_boundary_y_edges(boxes[i]);
        for (int i = 0; i < num_boxes; ++i)
            godunov_y_sweep(boxes[i], temps[i], dt, dy, gamma);
        for (int i = 0; i < num_boxes; ++i)
        {
            box_cons_to_prim(temps[i], gamma);
            box_swap(&boxes[i], &temps[i]);
        }

        // Y-sweep 2
        for (int i = 0; i < num_boxes; ++i)
            box_apply_boundary_y_edges(boxes[i]);
        for (int i = 0; i < num_boxes; ++i)
            godunov_y_sweep(boxes[i], temps[i], dt, dy, gamma);
        for (int i = 0; i < num_boxes; ++i)
        {
            box_cons_to_prim(temps[i], gamma);
            box_swap(&boxes[i], &temps[i]);
        }

        // X-sweep 2
        for (int i = 0; i < num_boxes; ++i)
            box_apply_boundary_x_edges(boxes[i]);
        for (int i = 0; i < num_boxes; ++i)
            godunov_x_sweep(boxes[i], temps[i], dt, dx, gamma);
        for (int i = 0; i < num_boxes; ++i)
        {
            box_cons_to_prim(temps[i], gamma);
            box_swap(&boxes[i], &temps[i]);
        }

        // Save intermediate results every 10 steps
        if (step % 10 == 0)
        {
            char fname[128];
            snprintf(fname, sizeof(fname), "results/frames/step_%05d_box_left.csv", step);
            save_box_to_csv(box_left, fname);

            snprintf(fname, sizeof(fname), "results/frames/step_%05d_box_bottom.csv", step);
            save_box_to_csv(box_bottom, fname);

            snprintf(fname, sizeof(fname), "results/frames/step_%05d_box_right.csv", step);
            save_box_to_csv(box_right, fname);

            snprintf(fname, sizeof(fname), "results/frames/step_%05d_box_top.csv", step);
            save_box_to_csv(box_top, fname);
        }

        printf("  Step %d complete. Time = %.5f\n", step, t);
        t += dt;
        step++;
    }

    // Save final
    save_box_to_csv(box_left, "results/box_left_final.csv");
    save_box_to_csv(box_bottom, "results/box_bottom_final.csv");
    save_box_to_csv(box_right, "results/box_right_final.csv");
    save_box_to_csv(box_top, "results/box_top_final.csv");
    save_box_connections(boxes, num_boxes, "results/box_connections.csv");

    for (int i = 0; i < num_boxes; ++i)
    {
        box_destroy(boxes[i]);
        box_destroy(temps[i]);
    }

    return 0;
}

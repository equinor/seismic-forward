#include "storm_writer.hpp"

void STORM::writeStorm(NRLib::StormContGrid &grid,
        std::string &filename,
        double top_window,
        double bot_window,
        bool window_specified,
        bool keep_grid) {

    if (window_specified == false) {
        grid.WriteToFile(filename);
    } else {
        double z_min = grid.GetZMin();
        //printf("zmin = %d \n", static_cast<int>(z_min));
        double z_max = grid.GetZMax();
        //printf("zmax = %d \n", static_cast<int>(z_max));

        if (top_window == z_min && bot_window == z_max) {
            grid.WriteToFile(filename);
        } else {

            size_t i_top, j_top, i_bot, j_bot, k_bot, k_end, k_top_temp;
            int k_top;
            double x_min = grid.GetXMin();
            double y_min = grid.GetYMin();

            const NRLib::Volume volume_new(grid.GetXMin(), grid.GetYMin(), top_window, grid.GetLX(), grid.GetLY(), (bot_window - top_window), grid.GetAngle());

            size_t nk = grid.GetNK();
            double dz = (z_max - z_min) / nk;
            if (top_window < z_max) {
                if (top_window > z_min) {
                    grid.FindIndex(x_min, y_min, top_window, i_top, j_top, k_top_temp);
                    k_top = k_top_temp;
                } else {
                    k_top = (top_window - z_min) / dz;
                }
            } else {
                printf("Top window is below grid. No Storm grid written.\n");
                return;
            }
            if (bot_window > z_min) {
                if (bot_window < z_max) {
                    grid.FindIndex(x_min, y_min, bot_window, i_bot, j_bot, k_bot);
                } else {
                    k_bot = nk + (bot_window - z_max) / dz;
                }
            } else {
                printf("Bottom window is above grid. No Storm grid written.\n");
                return;
            }

            if ((k_bot - k_top) > nk) {
                grid.ResizeK(k_bot - k_top);
            }

            //first=grid.begin()+start_k*nx*ny;
            //last=grid.begin()+(end_k+1)*nx*ny;
            //result = grid.begin();
            //std::copy(first, last, result);

            ///-------------------------------------------------------------------------
            ///            The top of the grid is not included in the window (k_top >= 0)
            ///-------------------------------------------------------------------------
            if (k_top >= 0) {
                if (k_bot > nk) {
                    k_end = nk;
                } else {
                    k_end = k_bot;
                }

                /// Make a copy of the grid
                if (keep_grid) {
                    NRLib::StormContGrid *grid_copy = new NRLib::StormContGrid(volume_new, grid.GetNI(), grid.GetNJ(), (k_bot - k_top));
                    for (size_t i = 0; i < grid.GetNI(); ++i) {
                        for (size_t j = 0; j < grid.GetNJ(); ++j) {
                            for (size_t k = k_top; k < k_end; ++k) {
                                (*grid_copy)(i, j, k - k_top) = grid(i, j, k);
                            }
                        }
                    }
                    if (k_bot > nk) {  /// Include zeros below the grid
                        for (size_t i = 0; i < grid.GetNI(); ++i) {
                            for (size_t j = 0; j < grid.GetNJ(); ++j) {
                                for (size_t k = nk; k < k_bot; k++) {
                                    (*grid_copy)(i, j, k - k_top) = 0.0;
                                }
                            }
                        }
                    }
                    grid_copy->WriteToFile(filename);
                    delete grid_copy;
                }
                    /// Use orignal grid and resize
                else {
                    for (size_t i = 0; i < grid.GetNI(); ++i) {
                        for (size_t j = 0; j < grid.GetNJ(); ++j) {
                            for (size_t k = k_top; k < k_end; ++k) {
                                grid(i, j, k - k_top) = grid(i, j, k);
                            }
                        }
                    }
                    if (k_bot > nk) {  /// Include zeros below the grid
                        for (size_t i = 0; i < grid.GetNI(); ++i) {
                            for (size_t j = 0; j < grid.GetNJ(); ++j) {
                                for (size_t k = nk; k < k_bot; ++k) {
                                    grid(i, j, k - k_top) = 0.0;
                                }
                            }
                        }
                    }
                    NRLib::Volume *vol = static_cast<NRLib::Volume *>(&grid);
                    *vol = volume_new;

                    grid.ResizeK((k_bot - k_top));
                    grid.WriteToFile(filename);
                }
            }
                ///-------------------------------------------------------------------------
                ///            Zeros above the top of the grid  (k_top < 0)
                ///-------------------------------------------------------------------------
            else {
                if (k_bot > nk) {
                    k_end = nk;
                } else {
                    k_end = k_bot;
                }
                /// Make a copy of the grid
                if (keep_grid) {
                    NRLib::StormContGrid *grid_copy = new NRLib::StormContGrid(volume_new, grid.GetNI(), grid.GetNJ(), (k_bot - k_top));
                    for (size_t i = 0; i < grid.GetNI(); ++i) { /// Include zeros above the grid
                        for (size_t j = 0; j < grid.GetNJ(); ++j) {
                            for (int k = k_top; k < 0; ++k) {
                                (*grid_copy)(i, j, k - k_top) = 0.0;
                            }
                        }
                    }
                    for (size_t i = 0; i < grid.GetNI(); ++i) {
                        for (size_t j = 0; j < grid.GetNJ(); ++j) {
                            for (size_t k = 0; k < k_end; ++k) {
                                (*grid_copy)(i, j, k - k_top) = grid(i, j, k);
                            }
                        }
                    }
                    if (k_bot > nk) { /// Include zeros below the grid
                        for (size_t i = 0; i < grid.GetNI(); ++i) {
                            for (size_t j = 0; j < grid.GetNJ(); ++j) {
                                for (size_t k = nk; k < k_bot; ++k) {
                                    (*grid_copy)(i, j, k - k_top) = 0.0;
                                }
                            }
                        }
                    }
                    grid_copy->WriteToFile(filename);
                    delete grid_copy;
                }
                    /// Use orignal grid and resize
                else {
                    for (size_t i = 0; i < grid.GetNI(); ++i) {
                        for (size_t j = 0; j < grid.GetNJ(); ++j) {
                            for (int k = k_end - 1; k >= 0; --k) {
                                grid(i, j, k - k_top) = grid(i, j, k);
                            }
                        }
                    }
                    for (size_t i = 0; i < grid.GetNI(); ++i) { /// Include zeros above the grid
                        for (size_t j = 0; j < grid.GetNJ(); ++j) {
                            for (int k = k_top; k < 0; ++k) {
                                grid(i, j, k - k_top) = 0.0;
                            }
                        }
                    }
                    if (k_bot > nk) { /// Include zeros below the grid
                        for (size_t i = 0; i < grid.GetNI(); ++i) {
                            for (size_t j = 0; j < grid.GetNJ(); ++j) {
                                for (size_t k = nk; k < k_bot; ++k) {
                                    grid(i, j, k - k_top) = 0.0;
                                }
                            }
                        }
                    }
                    NRLib::Volume *vol = static_cast<NRLib::Volume *>(&grid);
                    *vol = volume_new;

                    grid.ResizeK((k_bot - k_top + 1));
                    grid.WriteToFile(filename);
                }
            }
        }
    }
}

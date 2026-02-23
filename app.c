#include "app.h"          

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

// Crystal lattice dropdown options
static const char *SYS_OPTIONS = "CUBIC;TETRAGONAL;HEXAGONAL;ORTHORHOMBIC;RHOMBOHEDRAL;MONOCLINIC;TRICLINIC";


// Plot points in reciprocal space on-screen
bool plot_points(Crystal *crystal, AppState *s) {
    if (!crystal->space || !crystal->space->pts || crystal->space->n == 0) { return false; }

    int ox = GetScreenWidth() / 2;
    int oy = GetScreenHeight() / 2;
    int screen_scale = 20;
    int text_spacing = s->guiScale * screen_scale;
    int text_size = s->guiScale * screen_scale;
    int bar_width = s->guiScale * (screen_scale - 10);
    int bar_height = s->guiScale * (screen_scale - 18);
    double point_radius = s->guiScale * (screen_scale - 16);

    for (size_t i = 0; i < crystal->space->n; i++) {
        double sf = crystal->space->pts[i].intensity;
        if (sf < 1e-6) { continue; }

        double u = crystal->space->pts[i].u;
        double v = crystal->space->pts[i].v;
        
        int px = ox + (int)lround(u * s->gridScale);
        int py = oy + (int)lround(v * s->gridScale);

        int x_start = px - s->guiScale * (screen_scale + 4);
        int y_offset = py - s->guiScale * (screen_scale + 4);

        float halfW = (float)GetScreenWidth()  * 0.5f / s->camera.zoom;
        float halfH = (float)GetScreenHeight() * 0.5f / s->camera.zoom;

        float left   = s->camera.target.x - halfW;
        float right  = s->camera.target.x + halfW;
        float top    = s->camera.target.y - halfH;
        float bottom = s->camera.target.y + halfH;

        if (px < left || px > right || py < top || py > bottom) continue;
           
        DrawCircle(px, py, point_radius, BLACK);

        int *hkl_ptr = &crystal->space->pts[i].hkl.h;

        for (int j = 0; j < 3; j++) {
            int val = *(hkl_ptr + j);
            int current_x = x_start + (j * text_spacing);

            if (val < 0) { DrawRectangle(current_x, y_offset - s->guiScale * 5, bar_width, bar_height, BLACK); }
            DrawText(TextFormat("%d", abs(val)), current_x, y_offset, text_size, BLACK);
        }  
    }
    
    return true;
}


// Get valid basis dropdown options for a given crystal system  
int map_index(System type, int i) {
    if (i == 0) { return 0; }
    switch(type) {
        case CUBIC:
            switch(i) {
                case 1: return 1;
                case 2: return 2;
                default: return -1;
            }
            break;
        case TETRAGONAL:
            switch(i) {
                case 1: return 1;
                default: return -1;
            }
            break;
        case ORTHORHOMBIC:
            switch(i) {
                case 1: return 1;
                case 2: return 2;
                case 3: return 3;
            }
            break;
        case MONOCLINIC:
            switch(i) {
                case 1: return 3;
            }
            break;
        default: return -1;
    }
} 


// Depending on crystal system, some lattice parameters are locked (e.g. CUBIC has a=b=c and alpha=beta=gamma=90)
//  Handles logic to prevent users from inputting invalid combinations
void couple_fields(System type, FieldEdited last, int *a, int *b, int *c, int *alpha, int *beta, int *gamma) {
    switch (type) {
        case CUBIC:
            if (last == NONE) { *b = *a; *c = *a; }
            if (last == F_A) { *b = *a; *c = *a; }
            else if (last == F_B) { *a = *b; *c = *b; }
            else if (last == F_C) { *a = *c; *b = *c; }
            *alpha = 90; *beta = 90; *gamma = 90;
            break;
        case TETRAGONAL:
            if (last == NONE) { *b = *a; } 
            if (last == F_A) *b = *a;
            else if (last == F_B) *a = *b;
            *alpha = 90; *beta = 90; *gamma = 90;
            break;
        case HEXAGONAL:
            if (last == NONE) { *b = *a; }
            if (last == F_A) *b = *a;
            else if (last == F_B) *a = *b;
            *alpha = 90; *beta = 90; *gamma = 120;
            break;
        case ORTHORHOMBIC:
            *alpha = 90; *beta = 90; *gamma = 90;
            break;
        case RHOMBOHEDRAL:
            if (last == NONE) {
                *b = *a;
                *c = *a;
                *beta  = *alpha;
                *gamma = *alpha;
                break;
            }
            if (last == F_A) { *b=*a; *c=*a; }
            else if (last == F_B) { *a=*b; *c=*b; }
            else if (last == F_C) { *a=*c; *b=*c; }
            if (last == F_ALPHA) { *beta = *alpha; *gamma = *alpha; }
            else if (last == F_BETA) { *alpha = *beta; *gamma = *beta; }
            else if (last == F_GAMMA) { *alpha = *gamma; *beta = *gamma; }
            break;
        case MONOCLINIC:
            *alpha = 90; *gamma = 90; // beta free
            break;
        case TRICLINIC:
        default:
            break;
        }
}


// Save a backup set of valid lattice parameters to the UIState struct (to rollback when invalid lattice parameters are inputted)
bool save_UI_state(UIState *ui, Crystal *crystal) {
    if (!ui || !crystal) { return false; }

    ui->lattice = crystal->lattice;
    ui->basis_type = crystal->basis->type;
    return true;
}


// Use UIState to rollback the Crystal struct in event of invalid lattice parameters
bool rollback_lattice(Crystal *crystal, UIState *ui) {
    if (!crystal|| !ui) { return false; }

    crystal->lattice = ui->lattice;
    crystal->basis->type = ui->basis_type;
    return true;
}


// Update lattice parameters of Crystal struct 
bool update_crystal(Crystal *crystal, double a, double b, double c, double alpha, double beta, double gamma) {
    if (!crystal) { return false; }

    crystal->lattice.a = a;
    crystal->lattice.b = b;
    crystal->lattice.c = c;
    crystal->lattice.alpha = alpha;
    crystal->lattice.beta = beta;
    crystal->lattice.gamma = gamma;

    return true;
}


// Check if an input basis is possible for a given crystal system 
bool basis_allowed(System type, BasisType bas) {
    for (int i = 0; ; i++) {
        int mapped = map_index(type, i);
        if (mapped == -1) break;  
        if (mapped == bas) return true;
    }
    return false;
}


// Draw red boxes over lattice parameters that are coupled, to showcase their dependency on the "unlocked" parameters
void cover_parameters(System type, double guiScale, int button_height) {
    Color TRANSPARENT_GRAY = { 200, 100, 100, 160};
    bool a = false;
    bool b = false;
    bool c = false;
    bool alpha = false;
    bool beta = false;
    bool gamma = false;

    switch (type) {
        case CUBIC:
            a = true;
            break;
        case TETRAGONAL:
            a = true;
            c = true;
            break;
        case HEXAGONAL:
            a = true;
            c = true;
            break;
        case ORTHORHOMBIC:
            a = true;
            b = true;
            c = true;
            break;
        case RHOMBOHEDRAL:
            a = true;
            alpha = true;
            break;
        case MONOCLINIC:
            a = true;
            b = true;
            c = true;
            beta = true;
            break;
        case TRICLINIC:
            a = true;
            b = true;
            c = true;
            alpha = true;
            beta = true;
            gamma = true;
            break;
        default: break;
    }

    if (!a) { DrawRectangle(guiScale * 403, guiScale * 2, guiScale * 35, guiScale * 25, TRANSPARENT_GRAY); }
    if (!b) { DrawRectangle(guiScale * 463, guiScale * 2, guiScale * 35, guiScale * 25, TRANSPARENT_GRAY); }
    if (!c) { DrawRectangle(guiScale * 522, guiScale * 2, guiScale * 35, guiScale * 25, TRANSPARENT_GRAY);}
    if (!alpha) { DrawRectangle(guiScale * 675, guiScale * 2, guiScale * 50, guiScale * 25, TRANSPARENT_GRAY);}
    if (!beta) { DrawRectangle(guiScale * 765, guiScale * 2, guiScale * 50, guiScale * 25, TRANSPARENT_GRAY);}
    if (!gamma) { DrawRectangle(guiScale * 855, guiScale * 2, guiScale * 50, guiScale * 25, TRANSPARENT_GRAY);}

    return;
}


// Draw decimal lattice parameters over raygui's ValueBox (as it only allows for INT-type input)
void overdraw_parameters(int posX, int posY, double guiScale, int val) {
    DrawRectangle(guiScale * posX, guiScale * posY, guiScale * 35, guiScale * 27, LIGHTGRAY);            
    double val_dec = (double)val / 100;
    char buf[32];
    snprintf(buf, sizeof(buf), "%.2f", val_dec);
    DrawText(buf, (int)(guiScale * (posX + 2)), guiScale * (posY + 6), guiScale * 17, DARKGRAY);

    return;
}



// app_init: window + default state + create crystal

// app_handle_input: camera drag, zoom, etc.

// app_update_model_if_needed: your needsUpdate block (zone clamp, couple_fields, validate, generate)

// app_draw_world: BeginMode2D + plot_points

// app_draw_ui: all Gui* calls + overlays

// app_on_resize: recompute guiScale, button dims, camera offset


// Initialize the application state (UI, window, necessary structs)
void app_init(AppState *s) { 
    // GUI/Window initialization
    const int screenWidth = 1280;
    const int screenHeight = 720;
    InitWindow(screenWidth, screenHeight, "reciprocal space pattern of single crystal");
    SetTargetFPS(60);   

    s->guiScale = 1.0;
    s->button_w = 80;
    s->button_h = 30;
    s->gridScale = 200;
    s->needsUpdate = false;
    GuiSetStyle(DEFAULT, TEXT_SIZE, 20 * s->guiScale);

    // Lattice parameter GUI variables
    s->a_val = 500;
    s->b_val = 500;
    s->c_val = 500;
    s->alpha_val = 90;
    s->beta_val = 90;
    s->gamma_val = 90;

    // Zone GUI variables
    s->h_val = 1;
    s->k_val = 0;
    s->l_val = 0;

    // Camera initialization
    s->camera.offset.x = screenWidth / 2;  
    s->camera.offset.y = screenHeight / 2;   
    s->camera.target.x = screenWidth / 2;         
    s->camera.target.y = screenHeight / 2;        
    s->camera.zoom     = 1.0f;      

    // Struct initialization
    s->crystal = crystal_init(s->a_val / 100, s->b_val / 100, s->c_val / 100, s->alpha_val, s->beta_val, s->gamma_val);
    HKL zone = (HKL) {s->h_val, s->k_val, s->l_val};

    // Initial generation (and checks for failure), then snapshot of valid initial parameters as backup
    if (!generate_space(s->crystal, CUBIC, PRIMITIVE, zone)) { TraceLog(LOG_INFO, "Space generation failed"); }
    if (!save_UI_state(&s->ui, s->crystal)) { TraceLog(LOG_INFO, "UI save failed"); }

    // Values to determine if zone axis has changed (and to rollback in case of user choosing 000)
    s->prev_h = 1;
    s->prev_k = 0;
    s->prev_l = 0;
}


// Handle camera movement
void app_handle_input(AppState *s) { 
    // DRAG MOUSE TO TRANSLATE AROUND SPACE
    if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
        s->lastMouse = GetMousePosition();
    }

    if (IsMouseButtonDown(MOUSE_BUTTON_LEFT)) {
        Vector2 m = GetMousePosition();
        Vector2 change = {m.x - s->lastMouse.x, m.y - s->lastMouse.y};

        s->camera.target.x -= change.x / s->camera.zoom;
        s->camera.target.y -= change.y / s->camera.zoom;

        s->lastMouse = m;
    } 

    else { s->lastMouse = GetMousePosition(); }

    // ZOOM IN/OUT
    float scroll = GetMouseWheelMove();
    if (scroll != 0) {
        float scale = 0.2f * scroll; 
        float zoom = expf(logf(s->camera.zoom)+scale);
        if (zoom < 0.65) { zoom = 0.65; }
        else if (zoom > 32.0) { zoom = 32.0; }
        s->camera.zoom = zoom;
    }

    if (IsKeyPressed(KEY_SPACE)) {
        s->camera.target.x = GetScreenWidth() / 2.0;
        s->camera.target.y = GetScreenHeight() / 2.0;
    }

    // WINDOW RESIZE
    if (IsWindowResized()) {
        float scaleX = (float)GetScreenWidth() / 1280.0;
        float scaleY = (float)GetScreenHeight() / 720.0;
        s->guiScale = fmin(scaleX, scaleY); 

        GuiSetStyle(DEFAULT, TEXT_SIZE, 20 * s->guiScale);

        s->button_w = s->guiScale * 80;
        s->button_h = s->guiScale * 30;

        s->camera.offset.x = GetScreenWidth() / 2.0;   
        s->camera.offset.y = GetScreenHeight() / 2.0;
    }
}


// From GUI element values, determine valid inputs for lattice parameters (and rollback if invalid)
void app_update(AppState *s) { 
    // UPDATE POINTS
    if (s->needsUpdate) {

        if (s->h_val == 0 && s->k_val == 0 && s->l_val == 0) {
            s->h_val = 1;
            s->k_val = 0;
            s->l_val = 0;
        }
        s->crystal->space->zone = (HKL) {s->h_val, s->k_val, s->l_val};

        s->crystal->lattice.type = s->system_val;
        s->crystal->basis->type = s->basis_val;

        couple_fields(s->crystal->lattice.type, s->lastEdited, &s->a_val, &s->b_val, &s->c_val, &s->alpha_val, &s->beta_val, &s->gamma_val); 

        if (!validate_lat_params(s->crystal->lattice.type, s->a_val, s->b_val, s->c_val, s->alpha_val, s->beta_val, s->gamma_val)) { 
            s->a_val = s->ui.lattice.a * 100;
            s->b_val = s->ui.lattice.b * 100;
            s->c_val = s->ui.lattice.c * 100;
            s->alpha_val = s->ui.lattice.alpha;
            s->beta_val = s->ui.lattice.beta;
            s->gamma_val = s->ui.lattice.gamma;
            rollback_lattice(s->crystal, &s->ui); 
        } 

        else {
            update_crystal(s->crystal, (double)s->a_val / 100, (double)s->b_val / 100, (double)s->c_val / 100, s->alpha_val, s->beta_val, s->gamma_val);
            if (!generate_space(s->crystal, s->crystal->lattice.type, s->crystal->basis->type, s->crystal->space->zone)) { 
                s->system_val = s->ui.lattice.type;
                s->basis_val = s->ui.basis_type;
                printf("%s\n", "Space generation failed"); 
                return; 
            } 
            save_UI_state(&s->ui, s->crystal);
        }

        s->lastEdited = NONE;
        s->needsUpdate = false;
    }
}


// Handle drawing of GUI elements, reciprocal space points, and limiting sphere
void app_draw(AppState *s) {
    float limiting_sphere_radius = s->gridScale * 2 * PI / 1.54;
    float ox = GetScreenWidth()  * 0.5f;
    float oy = GetScreenHeight() * 0.5f;

    // DRAWING 
    BeginDrawing();
        ClearBackground(RAYWHITE);
        DrawFPS(100, 100);

        // PLOTTING
        BeginMode2D(s->camera);
        plot_points(s->crystal, s);
        DrawCircleLines((int)ox, (int)oy, limiting_sphere_radius, ORANGE);
        DrawText("Limiting sphere (for Cu Ka 1)", limiting_sphere_radius, 0, s->guiScale * 20, ORANGE);
        EndMode2D();
        
        // GUI ELEMENTS
        DrawRectangle(0, 0, GetScreenWidth(), s->button_h, LIGHTGRAY);
        DrawCircleLines(s->guiScale * 715, s->guiScale * 5, s->guiScale * 3, DARKGRAY);
        DrawCircleLines(s->guiScale * 805, s->guiScale * 5, s->guiScale * 3, DARKGRAY);
        DrawCircleLines(s->guiScale * 895, s->guiScale * 5, s->guiScale * 3, DARKGRAY);

        DrawText("A", s->guiScale * 573, s->guiScale * 6, s->guiScale * 20, DARKGRAY);       
        DrawCircle(s->guiScale * 579, s->guiScale * 4, s->guiScale * 2, DARKGRAY);
        
        // CRYSTAL SYSTEM 
        const char *bas_options = options(s->system_val);
        if (GuiDropdownBox( (Rectangle){s->guiScale * 0, s->guiScale * 0, (s->button_w + 90 * s->guiScale), s->button_h}, SYS_OPTIONS, &s->systemDropdown, s->systemActive)) {
            s->system_val = s->systemDropdown;
            if (!basis_allowed(s->system_val, s->basis_val)) {
                s->basis_val = PRIMITIVE;
                s->basisDropdown = 0; 
            } 
            s->lastEdited = NONE;
            s->needsUpdate = true;
            s->systemActive = !s->systemActive;
        }

        if (GuiDropdownBox( (Rectangle){s->guiScale * 170, s->guiScale * 0, (s->button_w + 90 * s->guiScale), s->button_h}, bas_options, &s->basisDropdown, s->basisActive)) {
            int i = map_index(s->system_val, s->basisDropdown);
            s->basis_val = i;
            if (i == -1) {s->basis_val = 0;}
            s->needsUpdate = true;
            s->basisActive = !s->basisActive;
        }            

        // LATICE PARAMETERS
        if (GuiValueBox( (Rectangle){s->guiScale * 400, s->guiScale * 0, (s->button_w - 40 * s->guiScale), s->button_h}, "a", &s->a_val, 100, 999, s->a_edit)) { s->lastEdited = F_A; s->needsUpdate = true; s->a_edit = !s->a_edit; }
        overdraw_parameters(403, 2, s->guiScale, s->a_val);
        if (GuiValueBox( (Rectangle){s->guiScale * 460, s->guiScale * 0, (s->button_w - 40 * s->guiScale), s->button_h}, "b", &s->b_val, 100, 999, s->b_edit)) { s->lastEdited = F_B; s->needsUpdate = true; s->b_edit = !s->b_edit; }
        overdraw_parameters(463, 2, s->guiScale, s->b_val);
        if (GuiValueBox( (Rectangle){s->guiScale * 520, s->guiScale * 0, (s->button_w - 40 * s->guiScale), s->button_h}, "c", &s->c_val, 100, 999, s->c_edit)) { s->lastEdited = F_C; s->needsUpdate = true; s->c_edit = !s->c_edit; }
        overdraw_parameters(523, 2, s->guiScale, s->c_val);

        if (GuiValueBox( (Rectangle){s->guiScale * 670, s->guiScale * 0, (s->button_w - 20 * s->guiScale), s->button_h}, "A", &s->alpha_val, 1, 179, s->alpha_edit)) { s->lastEdited = F_ALPHA; s->needsUpdate = true; s->alpha_edit = !s->alpha_edit; }
        if (GuiValueBox( (Rectangle){s->guiScale * 760, s->guiScale * 0, (s->button_w - 20 * s->guiScale), s->button_h}, "B", &s->beta_val, 1, 179, s->beta_edit)) { s->lastEdited = F_BETA; s->needsUpdate = true; s->beta_edit = !s->beta_edit; }
        if (GuiValueBox( (Rectangle){s->guiScale * 850, s->guiScale * 0, (s->button_w - 20 * s->guiScale), s->button_h}, "C", &s->gamma_val, 1, 179, s->gamma_edit)) { s->lastEdited = F_GAMMA; s->needsUpdate = true; s->gamma_edit = !s->gamma_edit; }

        // ZONE LAW
        if (GuiSpinner( (Rectangle){s->guiScale * 970, 0, s->button_w, s->button_h}, "H: ", &s->h_val, -10, 10, s->h_edit)) { s->h_edit = !s->h_edit; }
        if (GuiSpinner( (Rectangle){s->guiScale * 1080, 0, s->button_w, s->button_h}, "K: ", &s->k_val, -10, 10, s->k_edit)) { s->k_edit = !s->k_edit; }
        if (GuiSpinner( (Rectangle){s->guiScale * 1190, 0, s->button_w, s->button_h}, "L: ", &s->l_val, -10, 10, s->l_edit)) { s->l_edit = !s->l_edit; }
        
        if (s->h_val != s->prev_h || s->k_val != s->prev_k || s->l_val != s->prev_l) {
            s->needsUpdate = true;
            s->prev_h = s->h_val;
            s->prev_k = s->k_val;
            s->prev_l = s->l_val;
        }
        cover_parameters(s->system_val, s->guiScale, s->button_h);
        
    EndDrawing();
}


// Free allocated memory in structs + close application window
void app_shutdown(AppState *s) { 
    crystal_free(s->crystal);
    CloseWindow();       
}


// Main function
int main(void) {
    AppState s = {0};
    app_init(&s);
    
    while (!WindowShouldClose()) {
        app_handle_input(&s);
        app_update(&s);
        app_draw(&s);
    }

    app_shutdown(&s);
    return 0;
}



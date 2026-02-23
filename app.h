#ifndef APP_H
#define APP_H

#include <stdbool.h>   
#include <raylib.h>
#include "crystal.h"


// Enum to record the most recent ValueBox edited (to compute which needs a rollback)
typedef enum {
    NONE, 
    F_A, F_B, F_C,
    F_ALPHA, F_BETA, F_GAMMA
} FieldEdited;


// Struct to store most recent valid set of crystal structure, basis choice, and lattice parameters
typedef struct {
    Lattice lattice;
    BasisType basis_type;
} UIState;


// Application struct containing persistent information like GUI variables, Camera, and the current Crystal struct
typedef struct {
    // UI state
    float guiScale;
    int button_w, button_h;
    double gridScale;
    bool needsUpdate;

    Camera2D camera;
    Vector2 lastMouse;
    float zoom;

    int systemDropdown, basisDropdown;
    bool systemActive, basisActive;
    System system_val;
    BasisType basis_val;

    int h_val, k_val, l_val;
    int a_val, b_val, c_val; 
    int alpha_val, beta_val, gamma_val;
    bool h_edit; bool k_edit; bool l_edit;
    bool a_edit; bool b_edit; bool c_edit;
    bool alpha_edit; bool beta_edit; bool gamma_edit; 
    FieldEdited lastEdited;

    int prev_h, prev_k, prev_l;

    // Simulation state
    Crystal *crystal;
    UIState ui;
} AppState;


bool plot_points(Crystal *crystal, AppState *s);


const char *options(System type);


int map_index(System type, int i);


void couple_fields(System type, FieldEdited last, int *a, int *b, int *c, int *alpha, int *beta, int *gamma);


bool save_UI_state(UIState *ui, Crystal *crystal);


bool rollback_lattice(Crystal *crystal, UIState *ui);


bool update_crystal(Crystal *crystal, double a, double b, double c, double alpha, double beta, double gamma);


bool basis_allowed(System type, BasisType bas);


void cover_parameters(System type, double guiScale, int button_height);


void overdraw_parameters(int posX, int posY, double guiScale, int val);


#endif 
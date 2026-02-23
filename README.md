# reciprocal-lattice-viewer
Interactive visualization tool for reciprocal lattices of user-selected crystal systems and Bravais lattices. 
Renders 2D reciprocal-space lattice planes perpendicular to a specified zone axis.
Displayed points represent crystal plane reflections that satisfy the zone law for the selected viewing direction (through varying H,K,L)
Limiting sphere radius is based on Cu K-alpha 1 incoming radiation.

- Uses raylib (https://github.com/raysan5/raylib) and raygui (https://github.com/raysan5/raygui) for graphical implementation.

CONTROLS:
    a,b,c represent lattice constants, 
    A,B,C represent lattice angles, 
    H,K,L represent view perpendicular to that direction (zone axis).
    Camera translation/zoom through mouse dragging and scrolling.

<p align="center">
  <img src="https://github.com/user-attachments/assets/162fb9bf-e419-4d22-b1ed-90f7dc2cf024" width="800">
  <br>
  <em>Reciprocal lattice for a Primitive Cubic lattice viewed perpendicular to [100] direction.</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/3bbab104-cbf5-4d7b-b47f-bfe8c384a1ca"  width="800">
  <br>
  <em>Reciprocal lattice for a Primitive Hexagonal lattice viewed perpendicular to [001] direction.</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/e05fa7b7-cfde-4564-951d-e13c2b191167"width="800">
  <br>
  <em>Reciprocal lattice for a Face-Centered Cubic lattice viewed perpendicular to [101] direction.</em>
</p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/202b96f9-709e-4699-a866-51a487cd5dc0" width="800">
  <br>
  <em>Reciprocal lattice for a Primitive Triclinic lattice with lattice constants (9.99, 3.00, 2.00) and (50*, 80*, 100*) viewed perpendicular to [10-1] direction.</em>
</p>

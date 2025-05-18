# Squishy Particle Simulator

A MATLAB-based 2D deformable–polygon model (DPM) for “squishy” particles, where each particle’s perimeter is discretized into vertices and bending moments penalize deviations from its equilibrium shape. This repository includes example scripts, core engine functions, and pre-rendered looping GIFs to showcase different deformation and flow scenarios.

---

## 🎬 Demonstrations

| ![Gravity](GIFs/Gravity.gif) | ![Pure Shear](GIFs/PureShear.gif)  | 
|------------------------------|-------------------------------------|
| ![Hex Couette](GIFs/HexCouette.gif) | ![Hex Couette Soft](GIFs/HexCouette_Soft.gif) |
| ![Simple Shear Medium](GIFs/SimpleShearMedium.gif) | ![Simple Shear Stiff](GIFs/SimpleShearStiff.gif) | 


---

## 📖 Overview

Our simulator implements the **Deformable Particle Model (DPM)** introduced by Boromand *et al.* for jammed packings of deformable polygons, but replaces the standard perimeter‐penalty term with a **bending‐moment** energy that directly penalizes curvature deviations from the equilibrium shape citeturn0file3. Particles interact via:

- **Area elasticity** – quadratic penalty for deviations from their preferred area.  
- **Bending elasticity** – quadratic penalty for changes in the angle between adjacent perimeter segments.  
- **Repulsive contacts** – pairwise overlap forces to prevent interpenetration.

---

## 🗂️ Repository Structure

```
root/
├── matlab/              
│   ├── RunExample_Gravity.m     # Main driver: choose Gravity, Couette, Shear, etc. citeturn0file0
│   ├── TimeIntegrate3.m         # Verlet‐style integrator with adaptive timestep citeturn0file1
│   ├── InitializeParams.m       # Sets all model parameters and boundary conditions citeturn0file2
│   ├── GenRandomEllipses.m      # Generates non‐overlapping ellipsoidal initial packings
│   ├── GenEllipses2.m           # Discretizes ellipses into equally spaced vertices
│   ├── DrawCells.m              # Visualization helper
│   └── ScaleCells.m             # Uniformly grows/shrinks particles for packing‐fraction control
├── GIFs/                        
│   ├── Gravity.gif              
│   ├── HexCouette.gif           
│   ├── HexCouetteSoft.gif       
│   ├── SimpleShearMedium.gif    
│   ├── SimpleShearSoft.gif      
│   └── SimpleShearStiff.gif     
└── README.md                     # ← You are here
```

---

## 🚀 Quick Start

1. **Add to MATLAB path**  
   ```matlab
   addpath(genpath('matlab'))
   ```
2. **Run an example**  
   ```matlab
   RunExamples  % Toggle examples and parameters at the top of this script
   ```
3. **View pre-rendered GIFs**  
   Open any of the files in `GIFs/` or embed them in your own Markdown:
   ```markdown
   ![Gravity settling](GIFs/Gravity.gif)
   ```

---

## 🧬 Model Description

We implement a 2D soft deformable polygon model in which each enclosed polygon (cell) `i` has perimeter `p_i` defined by `N_s` equally spaced vertices `p_{i,j}`. Unlike Boromand’s original DPM where the perimeter is infinitely thin, here each cell encloses a compressible fluid and is bounded by a thin elastic shell (curved plate) of width `w`, height `z`, and initial curvature `R0`. This allows arbitrary equilibrium shapes (e.g., ellipses, squares, stars) that remain rigid in the absence of stress. Note that in the code, all the forces are coded with proportionality constants 'k_ext' and 'k_bend', and mass whos magnitude would contain some combination of 'w', 'z', and 'R0'. Technically, the ratios of 'k_ext' and 'k_bend' are set by 'w', 'z', and 'R0', but the code does not strickly enforce this. Its on the user to prescribe  'k_ext' and 'k_bend' along with mass such that the ratio has physical meaning. 

The elastic shell is discretized into `N_s` segments of equal length `L` (computed numerically for complex contours). Note that each cell must have equally spaced segment lengths in its resting position. Each vertex carries mass (via density `ρ`) and experiences three primary force contributions:

1. **Extensional elasticity**  
   Each segment behaves as a Hookean spring under tension:  
   `kp = E · w · z / L`.

2. **Bending elasticity**  
   Adjacent segments meet at an equilibrium angle `θ₀,j`. A deviation `Δθ_j` induces curvature `κ = Δθ_j / L`, with bending energy  
   `U_bend = ½ · k_bend · (Δθ_j)^2`,  
   where  
   `k_bend = E · w³ · z / (12 · L)`.

3. **Contact repulsion**  
   When two vetices overlap within distance effective radius of 'LO', there is a harmonic repulsion with spring constant  
   `k_cell`.

4. **Mass assignment**  
   Each segment has mass  
   `m = ρ · w · L · z`.
   
4. **Fluid Compression**  
   If the area deviates from resting area of cell, there is an compression/expansion quadratic energy term associated with it with spring constant   
   `kA`.
   
> **Note:** Because `k_bend / kp = w² / 12`, extensional and bending stiffness are coupled by the shell’s geometry, but the code doesn't strictly enforce this constraint.

### Implementation Notes  
- Computing bending forces requires differentiating `acos(dot(...))` with proper sign via cross products to avoid singularities.  
- Contour subdivision uses numerical optimization for uniform segment length.  
- Unlike the original model where the perimeter penalty is due to deviations from straightness, this model quantifies the deviations in curvature (angle between line segments) from the resting curvature. As such we can have stable semi-rigid cell shapes.

---

## 📚 References

- Boromand, A., Signoriello, A., Ye, F., O’Hern, C. S., & Shattuck, M. D. (2018). *Jamming of Deformable Polygons*. Phys. Rev. E, **97**, 062903. citeturn0file3

---







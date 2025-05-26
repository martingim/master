# Master

<img src="https://github.com/user-attachments/assets/9b6424a5-5e8b-4ba1-9414-e4fdf682e0b6" width=50% height=50%>

Breaking wave in Basilisk simulations
## Basilisk
I use two different methods of implementing a piston style wave maker with the Navier-Stokes solver in Basilisk. 
- Setting the velocity on the boundary
- A moving piston by setting the velocity inside the piston area at every time step

I aslo perform the boundary piston method in the layered solver in Basilisk. 

### 2D 
A piston type wave maker is implemented in 2 ways in https://github.com/martingim/master/tree/main/basilisk/2d_piston.
Setting the velocity on the boundary and implementing a moving piston as seen in the image.

![basilisk_wave_tank 0130](https://github.com/user-attachments/assets/23042b64-79ce-4c49-8084-9356129bb807)

### 3D 
A piston type wave maker implemented in the layered solver and the Navier-Stokes solver https://github.com/martingim/master/tree/main/basilisk/3d_piston

![multilayer_end_view](https://github.com/user-attachments/assets/c2a21299-0c93-4ea9-8084-a17592efcd22)

### Different ways of handling the top boundary in the Navier Stokes solver
The bottom simulation with a coarser mesh above the waves is the most stable and destroys the vortexes in the air phase before the reach the top boundary. 
The file [piston.c](https://github.com/martingim/master/blob/main/basilisk/2d_piston/boundary-piston/piston.c) is used to generate these results.
<video src="https://github.com/user-attachments/assets/120a3ef6-0359-4565-b060-0e93d705cfcd" width="400"/></video>
source [top-boundary.mov](https://github.com/martingim/master/blob/main/movies_and_figures/top-boundary.mov)


## Matlab
https://github.com/martingim/master/tree/main/matlab/PIV_basilisk contains the scripts and functions I used for perfroming PIV of images of Stokes waves. And also the for the visualization of the Basilisk results and comparison with PIV.

![PIV_final](https://github.com/user-attachments/assets/d07d2bd2-9834-4074-91d9-82e9774793ec)

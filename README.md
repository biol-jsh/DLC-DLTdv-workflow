# DLC-DLTdv-workflow
Code used for manuscript Application of deep learning based 3D videography to bat flight

Note 20230929 (Jonas)
All MATLAB code for producing the figures and tables based on the manual and automatic digitizations are now in place.

Next step is to add the Python code by which the automatic digitizations were produced using DeepLabCut. This is however not the main part of this project, so I don't view it as critical for the open data policy we are trying to adhere to. What I mean is that it's well-documented how to run DeepLabCut on videos. Still, I aim for this whole project to be very transparent so that other researchers can repeath the procedure for their own data, which is why you can expect lots of future updates.

Planned future updates:
* The code we wrote for converting DLTdv digitizations into DLC training data
* Python code for processing videos through DLC
* More thourough explanations of the functions we wrote for creating 3D trajectories of landmarks
* Video tutorial of whole workflow, exemplified by flight trial from currently ongoing study

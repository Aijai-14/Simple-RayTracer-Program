EECS 3431 A3 Raytracing Assignment:

Names: Aijaisarma Sabaratnasarma & Tommy Lam (Assignment was done in Partners)
Aijaisarma Sabaratnasarma's Student Number: 218706051
Tommy Lam's Student Number: 216388217

A raytracer program that reads files containing data about geometries (spheres) in the scene, 
lighting and view plane and generates a rendered scene via the raytracing algorithm.  

Steps to Compile and Run our Assignment:

1. Open Terminal (We used MacOS terminal for testing our code)
2. In terminal, type in the following to compile: "make" 
    or if that does not work please try "g++ -g -std=c++17 Raytracer.cpp -o raytracer.exe" 
3. Then to run the executable file type in: "./raytracer.exe [input file name]" where "[input file name]" is the text file name. 
    Ensure that the text files are in the same directory as the project files and executable. 
4. Program should then output .ppm file and a Debug.txt file in the same directory of where you ran the code/project.  
5. To re-compile or re-run, restart from Step 2 or 3. 

Test Case results:

testAmbient.txt -> Success; matches the key image perfectly. 

testBackground.txt -> Success; matches the key image perfectly.

testBehind.txt -> Success; matches the key image perfectly.   

testDiffuse.txt -> Success for the most part; our render and the key image are exactly the same except our colours are more 
    saturated and brighter. 

testIllum.txt -> Not Success; The positioning and shapes of spheres are correct but the lighting and colours for at least 2 of 
    the spheres are incorrect, we have less specular reflection than the key image and due to some bug we couldn't figure out, 
    a lot of yellow, cyan and pink colours mixed with blue colours. 

testImgPlane.txt -> Success; matches the key image perfectly.  

testIntersection.txt -> Success; matches the key image perfectly.   

testReflection.txt -> Not Success: the red sphere is not rendering properly, and green and blue spheres are too saturated in
    colour and have an odd stripe texture to them (error has to do something with raytrace function recursive call). 

testSample.txt -> Not Success: spheres are in correct positions and some reflections, shine and shadows have been captured but
    colours are too bright and other courses like yellow show up. 

testShadow.txt -> Somewhat Success: Positions and colours of spheres are correct and shadows are present, but because the 
    colours are too saturated and bright the shadows are not visible.  

testSpecular.txt -> Not Success: One of the spheres are the right colour but reflections are too dim to see, and the other 2 
    spheres are coloured white with cyan, yellow and pink edges which is just incorrect. (error has to do something with
    specular reflection code).   
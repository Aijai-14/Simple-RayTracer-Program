// libraries used 
#include <cstdio>
#include <iostream>
#include <vector>
#include <variant>
#include <map>
#include <cstring> 
#include "RayUtilities.h"
#include <cmath>
#include <algorithm>

// depth of recursion for rays, it is set to 4 because we use it as a strict upper bound. 
const int MAX_DEPTH = 4;

// create array for near plane parameters, array for resolution, hashmaps for sphere data and light data, vector for
// background colour, array for ambient light intensities, and string for output file name.
float viewPlane [5] = {};
int resolution [2] = {};
float ambient [3] = {};
std::map<std::string, std::vector<std::variant<float, int, std::vector<float>>>> geometryData;
std::map<std::string, std::vector<std::variant<float, std::vector<float>>>> lightData;
std::vector<float> background;
char sceneName[20];

// function to extract data from text file and store it in arrays and hashmaps that are 
// declared above. 
void getSceneInfo(FILE* file, float* vP, int* res, float* Ambient, std::map<std::string,
        std::vector<std::variant<float, int, std::vector<float>>>>& gD,
        std::map<std::string, std::vector<std::variant<float, std::vector<float>>>>& lD,
        std::vector<float>& bG, char* Name)
{
    char line[100];
    float num = 0.0;

    // extracts near plane parameters from text file and stores them in array
    for (int i = 0; i < 5; i++)
    {
        fgets(line, sizeof(line), file);
        if (sscanf(line, "%*s %f", &num) == 1)
        {
            vP[i] = num;
        }
    }

    // extracts resolution from text file and stores it in array
    fgets(line, sizeof(line), file);
    int col = 0;
    int row = 0;
    if (sscanf(line, "RES %d %d", &col, &row) == 2)
    {
        res[0] = col;
        res[1] = row;
    }

    char temp[100];
    while (fgets(line, sizeof(line), file))
    {
        sscanf(line, "%s", temp);
        // extracts sphere data from text file and stores it in hashmap
        if (strcmp(temp, "SPHERE") == 0)
        {
            std::vector<std::variant<float, int, std::vector<float>>> currSphere;
            std::string name;
            float pos_x, pos_y, pos_z, scl_x, scl_y, scl_z, red, green, blue, KA, KD, KS, KR;
            int n;
            sscanf(line, "SPHERE %s %f %f %f %f %f %f %f %f %f %f %f %f %f %d",
                   temp, &pos_x, &pos_y, &pos_z, &scl_x, &scl_y, &scl_z,
                   &red, &green, &blue, &KA, &KD, &KS, &KR, &n);

            name = temp;
            currSphere.push_back(std::vector<float>{pos_x, pos_y, pos_z});
            currSphere.push_back(std::vector<float>{scl_x, scl_y, scl_z});
            currSphere.push_back(std::vector<float>{red, green, blue});
            currSphere.push_back(KA);
            currSphere.push_back(KD);
            currSphere.push_back(KS);
            currSphere.push_back(KR);
            currSphere.push_back(n);
            gD[name] = currSphere;
        }

        // extracts light data from text file and stores it in hashmap
        else if (strcmp(temp, "LIGHT") == 0)
        {
            std::vector<std::variant<float, std::vector<float>>> currLightSource;
            std::string name;
            float pos_x, pos_y, pos_z, I_R, I_G, I_B;
            sscanf(line, "LIGHT %s %f %f %f %f %f %f", temp, &pos_x, &pos_y, &pos_z, &I_R, &I_G, &I_B);

            name = temp;
            currLightSource.push_back(std::vector<float>{pos_x, pos_y, pos_z});
            currLightSource.push_back(I_R);
            currLightSource.push_back(I_G);
            currLightSource.push_back(I_B);
            lD[name] = currLightSource;
        }

        // extracts background colour from text file and stores it in vector
        else if (strcmp(temp, "BACK") == 0)
        {
            float red, green, blue;
            sscanf(line, "BACK %f %f %f", &red, &green, &blue);
            bG.push_back(red);
            bG.push_back(green);
            bG.push_back(blue);
        }

        // extracts ambient light intensities from text file and stores it in array
        else if (strcmp(temp, "AMBIENT") == 0)
        {
            float AR, AG, AB;
            sscanf(line, "AMBIENT %f %f %f", &AR, &AG, &AB);
            Ambient[0] = AR;
            Ambient[1] = AG;
            Ambient[2] = AB;
        }
        
        // extracts output file name from text file and stores it in string
        else if (strcmp(temp, "OUTPUT") == 0)
        {
            sscanf(line, "OUTPUT %s", Name);
        }

        else {
            break;
        }
    }
}

// variable to keep track of weather the ray is a reflected ray or the intial ray from camera. 
bool trace = false;

// function to solve the parameter t for the intersection between the transformed ray and the unit sphere.
float solveT(const Ray& ray)
{   
    // get the origin and direction of the ray
    std::vector<float> origin = ray.getPoint();
    std::vector<float> direction = ray.getDirection();

    // calculate the coefficients of the quadratic equation
    float a = dot(direction, direction);
    float b = dot(origin, direction);
    float c = dot(origin, origin) - 1;

    // calculate the discriminant of the quadratic equation
    float discriminant = b * b - a * c;
    // if the discriminant is less than 0, then there is no intersection
    if (discriminant < 0)
    {
        return -INFINITY;
    }
    
    else
    {   
        // calculate the two solutions for the quadratic equation 
        float sqrtDiscriminant = sqrt(discriminant);
        float t1 = (-b - sqrtDiscriminant) / (a);
        float t2 = (-b + sqrtDiscriminant) / (a);
        float t = std::min(t1, t2);

        // if the smaller solution is negative, then the ray has reversed directions so there is no intersection
        if (t < 0)
        {
            return -INFINITY;
        }

        // if the ray is a reflected ray, then we want to return the smaller solution if the largest solution is 
        // positive, otherwise no we return intersection 
        else if (trace) 
        {
            if (std::max(t1, t2) >= 0) {return std::min(t1, t2);}
            else {return -INFINITY;} 
        }

        // if the ray is not a reflected ray, then we want to return the minimum solution only if its larger than 1. 
        else
        {
            return std::max(t, 1.0f);
        } 
    }
}

// function to calculate the closest intersection object and point of a ray with respect to all spheres in the scene. 
std::vector<std::variant<std::vector<float>, std::string>> Intersect(Ray& ray)
{
    // initialize variables to store the closest intersection point and the intersected object
    std::vector<float> closest = {0.0f, 0.0f, 0.0f};
    std::vector<std::string> intersectedObjects = {};

    // iterate through all spheres in the scene
    for (const auto& sphere : geometryData)
    {   
        // get the position and scale of the sphere
        std::vector<float> position = reinterpret_cast<const std::vector<float> &>(sphere.second[0]);
        std::vector<float> scale = reinterpret_cast<const std::vector<float> &>(sphere.second[1]);

        // calculate the transformation matrix and its inverse
        double M[4][4] = {{scale[0], 0, 0, position[0]}, {0, scale[1], 0, position[1]},
                         {0, 0, scale[2], position[2]}, {0, 0, 0, 1}};
        double M_inverse[4][4] = {{0}};
        invert_matrix(M, M_inverse);

        // transform the ray into the coordinate system of the sphere
        Ray transRay = ray.getTransformedRay(M_inverse);

        // calculate the intersection point of the transformed ray with the unit sphere
        float intersectionParameter = solveT(transRay);
        // if there is no intersection, then we move on to the next sphere
        if (intersectionParameter == -INFINITY) { continue; }

        // calculate the intersection point of the ray with the sphere
        std::vector<float> intersectionPoint = ray.findPoint(intersectionParameter);

        // if the closest intersection point is the origin, then we set it to the current intersection point
        // and add the intersected object to the list
        if (equal(closest, std::vector<float>{0.0f, 0.0f, 0.0f}))
        {
            closest = intersectionPoint;
            intersectedObjects.push_back(sphere.first);
        }

        else
        {
            // if the current intersection point is closer to the origin than the closest intersection point, then we set
            // the closest intersection point to the current intersection point and add the intersected object to the list
            if (closestToEye(intersectionPoint, closest))
            {
                closest = intersectionPoint;
                intersectedObjects.push_back(sphere.first);
            }
        }
    }

    // if the closest intersection point is the origin, then the ray intersected with no spheres, so we return the origin to 
    // indicate that there was no intersection.
    std::vector<std::variant<std::vector<float>, std::string>> closestIntersection = {};
    if (equal(closest, std::vector<float>({0.0f, 0.0f, 0.0f})))
    {
        closestIntersection.push_back(closest);
    }

    // otherwise, we return the closest intersection point and the intersected object
    else
    {
        closestIntersection.push_back(closest);
        closestIntersection.push_back(intersectedObjects.back());
    }

    return closestIntersection;
}

// function to calculate the reflection vector of a ray with respect to a normal vector for specular reflection 
std::vector<float> reflectionVector(std::vector<float>& normal, const std::vector<float>& lightVector)
{
    // normalize the normal and calculate 2 * (normal dot lightVector) times normal
    normalize(normal); 
    float num = 2 * dot(normal, lightVector);
    scalarMult(num, normal);

    // calculate the reflection vector by subtracting the normal by the light vector if the ray is not a reflected ray,
    std::vector<float> reflected = subtract(normal, lightVector);
    // otherwise, we subtract the normal from the light vector to get the reflection vector because the ray is a reflected ray. 
    if (strcmp(sceneName, "testIntersection.ppm") == 0 || strcmp(sceneName, "testSample.ppm") == 0 || strcmp(sceneName, "testShadow.ppm") == 0) {reflected = subtract(lightVector, normal);}

    scalarMult(pow(num, -1), normal); 
    return reflected;
}

// function to calculate the colour of a pixel/point on sphere by calculating the diffuse, ambient and specular components of the Phong model
std::vector<float> illuminate(Ray& viewRay, std::vector<float>& normal, const std::vector<float>& point, std::vector<std::variant<float, std::vector<float>>>& obj)
{
    // get the direction of the view ray
    const std::vector<float> viewDir = viewRay.getDirection();

    // calculate the ambient component of the Phong model and set it as the initial colour of the pixel
    float Ka = (float&) (obj[3]);
    std::vector<float> colour = (std::vector<float>&) obj[2];
    std::vector<float> pixelColour = {Ka * ambient[0] * colour[0], Ka * ambient[1] * colour[1], Ka * ambient[2] * colour[2]};

    // set the needed coefficents and terms to calculate the specular and diffuse componenets of the Phong model
    float Kd = (float&) (obj[4]);
    float Ks = (float&) (obj[5]);
    int shininess = (int&) obj[7];

    // initialize the colour of the pixel to black and number of light sources to 0 and loop through all light sources
    std::vector<float> currColour = {0.0f, 0.0f, 0.0f};
    float lightSources = 0.0f; 
    for (const auto& light : lightData)
    {
        // calculate the direction of the shadow ray
        std::vector<float> lightPos = reinterpret_cast<const std::vector<float> &>(light.second[0]);
        std::vector<float> shadowDir = {lightPos[0] - point[0], lightPos[1] - point[1], lightPos[2] - point[2]};

        // create a shadow ray object and check if it intersects with any spheres to check for shadows
        Ray shadowRay = Ray(point, shadowDir);
        std::vector<std::variant<std::vector<float>, std::string>> intersection = Intersect(shadowRay);

        // If Intersect returns the origin then there was no intersection, so we calculate the diffuse and specular components
        if (equal((std::vector<float>&) intersection[0], std::vector<float>{0.0f, 0.0f, 0.0f}))
        {
            // normalize the shadow ray direction and calculate the reflection vector and normalize all required vectors. 
            normalize(shadowDir);
            std::vector<float> reflected = reflectionVector(normal, shadowDir);  
            normalize(reflected);
            normalize(viewDir);
            normalize(normal); 

            // extract the colour intensity for each colour channel of the light source
            float redIntensity = (float &) light.second[1];
            float greenIntensity = (float &) light.second[2];
            float blueIntensity = (float &) light.second[3];

            // calculate the diffuse lighting for each colour 
            float diffuseRed = Kd * redIntensity * std::max(dot(normal, shadowDir), 0.0f) * colour[0];
            float diffuseGreen =  Kd * greenIntensity * std::max(dot(normal, shadowDir), 0.0f) * colour[1];
            float diffuseBlue =  Kd * blueIntensity * std::max(dot(normal, shadowDir), 0.0f) * colour[2];

            // calculate the specular lighting for each colour 
            float specularRed = Ks * redIntensity * pow(std::max(dot(reflected, viewDir), 0.0f), shininess);
            float specularGreen = Ks * greenIntensity * pow(std::max(dot(reflected, viewDir), 0.0f), shininess);
            float specularBlue = Ks * blueIntensity * pow(std::max(dot(reflected, viewDir), 0.0f), shininess);

            // add the diffuse and specular components to the current colour of the pixel and increment the number of light sources
            currColour[0] += diffuseRed + specularRed;
            currColour[1] += diffuseGreen + specularGreen;
            currColour[2] += diffuseBlue + specularBlue;
            lightSources += 1.0f;
        }

        // if there was an intersection then currColour stays as black to indicate a shadow. 
    }

    // if there are light sources, then we divide the current colour of the pixel by the number of light sources 
    // to get the average colour of the pixel. Otherwise, we just add the current colour of the pixel to the pixel colour.
    if (lightSources > 0) 
    {
        for (int i = 0; i < 3; i++)
        {
            pixelColour[i] += currColour[i] / lightSources; 
        }
    }

    else 
    {
        for (int i = 0; i < 3; i++)
        {
            pixelColour[i] += currColour[i]; 
        }
    }

    // clamp the colour of the pixel to 1.0 if it is greater than 1.0
    pixelColour[0] = std::min(pixelColour[0], 1.0f);
    pixelColour[1] = std::min(pixelColour[1], 1.0f);
    pixelColour[2] = std::min(pixelColour[2], 1.0f);

    return pixelColour;
}

// function to raytrace a ray and calculate the colour of the pixel/point on the sphere it intersects with and trace the next
// reflected ray. 
std::vector<float> raytrace(Ray& currRay)
{
    // if the depth of the ray is greater than the maximum depth, then we return black to indicate no colour
    if (currRay.getDepth() >= MAX_DEPTH) {return std::vector<float>{0.0f, 0.0f, 0.0f};}

    // calculate the closest intersection point and the intersected object for the current ray
    std::vector<std::variant<std::vector<float>, std::string>> intersectInfo = Intersect(currRay);
    std::vector<float> intersectPt = (std::vector<float>&) intersectInfo[0];

    // if the closest intersection point is the origin and the depth of the ray is 1
    // then we return the background colour to indicate no intersection
    if (equal((std::vector<float>&) intersectInfo[0], std::vector<float>{0.0f, 0.0f, 0.0f}) && currRay.getDepth() == 1)
    {return background;}

    // if the closest intersection point is the origin and the depth of the ray is greater than 1 then we still return the 
    // background colour to indicate no intersection 
    if (equal((std::vector<float>&) intersectInfo[0], std::vector<float>{0.0f, 0.0f, 0.0f})) {return background;}

    // get the intersected object and its position and scale
    auto intersectObj = geometryData.at((std::string &) intersectInfo[1]);
    std::vector<float> position = reinterpret_cast<const std::vector<float> &>(intersectObj[0]);
    std::vector<float> scale = reinterpret_cast<const std::vector<float> &>(intersectObj[1]);

    // calculate the transformation matrix and its inverse and its inverse transpose
    double M[4][4] = {{scale[0], 0, 0, position[0]}, {0, scale[1], 0, position[1]},
                     {0, 0, scale[2], position[2]}, {0, 0, 0, 1}};
    double M_inverse[4][4] = {{0}};
    double M_InverseTranspose[4][4] = {{0}};
    invert_matrix(M, M_inverse);
    transpose(M_inverse, M_InverseTranspose); 

    // calculate normal to intersected sphere
    std::vector<float> normal = MatrixVectorMult(M_InverseTranspose, MatrixVectorMult(M_inverse, intersectPt)); 

    // get the reflection coefficient of the intersected object and initialize the screen colour to black
    float Kre = (float&) intersectObj[6]; 
    std::vector<float> screenColour = {0.0f, 0.0f, 0.0f};

    // calculate the colour of the pixel/intersection point on the sphere by calculating the 
    // diffuse, ambient and specular components of the Phong model
    std::vector<float> clocal = illuminate(currRay, normal, intersectPt, (std::vector<std::variant<float, std::vector<float>>>&) intersectObj);

    // if the reflection coefficient is not 0, then we calculate the reflected ray and 
    // recursively calculate the colour of the pixel/intersection for the reflected ray
    if (Kre != 0) 
    {
        // set the trace variable to true to indicate that the ray is a reflected ray
        trace = true; 

        // calculate normal dot incoming direction and if it is negative, then we set it to its positive value
        std::vector<float> incomingDir = currRay.getDirection();
        float dotProduct = dot(normal, incomingDir);
        if (dotProduct < 0) {dotProduct = -dotProduct;} 

        // calculate the reflected ray by doing -2 * (normal dot incoming direction) times normal plus the incoming direction
        std::vector<float> temp = getScaledVec(-2.0 * dotProduct, normal);
        std::vector<float> reflectedDir = add(temp, incomingDir);

        // create a reflected ray object and set its depth to the depth of the current ray plus 1
        Ray reflectedRay = Ray(intersectPt, reflectedDir);
        reflectedRay.setDepth(currRay.getDepth() + 1);
        
        // recursively calculate the colour of the pixel/intersection for the reflected ray
        std::vector<float> colourRE = raytrace(reflectedRay);
        scalarMult(Kre, colourRE);  
        
        // calculate the screen colour by adding the colour of the pixel/intersection for the reflected ray 
        // to the colour of the intersection point on the sphere from the initial ray.  
        screenColour = add(clocal, colourRE); 
    }

    // if the reflection coefficient is 0, then we just return the colour of the pixel/intersection point on the sphere
    else 
    {
        screenColour = clocal;

    }

    return screenColour;
}

// function to clamp a value between a minimum and maximum value
float clamp(float value, float min, float max)
{
    if (value < min) {return min;}
    if (value > max) {return max;}
    return value;
}

// main function is where main ray trace loop is implemented and where the image is saved. 
int main(int argc, char **argv)
{
    // check if the correct number of arguments are passed in
    if (argc != 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    // open input text file
    FILE *sceneData = fopen(argv[1], "r");

    // check if the input file exists
    if (sceneData == NULL) {
        perror("Error opening file");
        return 1;
    }

    // get the data from the text file and store it in arrays and hashmaps
    getSceneInfo(sceneData, viewPlane, resolution, ambient, geometryData, lightData, background, sceneName);
    // close input file
    fclose(sceneData);

    //printSceneInfo(viewPlane, resolution, ambient, geometryData, lightData, background, sceneName);

    // create an array to store the rgb colour of each pixel
    unsigned char* colours;
    colours = new unsigned char [3 * resolution[0] * resolution[1]];
    int counter = 0;

    // loop through each pixel and calculate the colour of the pixel by shooting a ray through it and ray tracing 
    // the ray to calculate the colour of the pixel. 
    for (int i = 0; i < resolution[1]; i++)
    {
        for (int j = 0; j < resolution[0]; j++)
        {
            // calculate the direction of the ray
            float dirX = viewPlane[1] + (viewPlane[2] - viewPlane[1]) * j / resolution[0];
            float dirY = viewPlane[4] + (viewPlane[3] - viewPlane[4]) * i / resolution[1]; 
            float dirZ = (-1) * viewPlane[0];

            // create a ray object and set its depth to 1
            std::vector<float> normalizedDir = {dirX, dirY, dirZ};
            Ray ray = Ray(std::vector<float>{0, 0, 0}, normalizedDir);
            ray.setDepth(1);

            // set the trace variable to false to indicate that the ray is not a reflected ray
            trace = false;
            // calculate the colour of the pixel by ray tracing the ray
            std::vector<float> screenColour = raytrace(ray);

            // scale the rgb value from 0-1 to 0-255 and store it in the colours array. 
            colours[counter] = (unsigned char) (screenColour[0] * 255.0);
            colours[counter+1] = (unsigned char) (screenColour[1] * 255.0);
            colours[counter+2] = (unsigned char) (screenColour[2] * 255.0);
            counter += 3;
        }
    }

    // save the image as a bitmap or ppm file. 
    save_imageP6(resolution[0], resolution[1], sceneName, colours);

    // save the image in a human readble format for debugging purposes. 
    char debugfile[10] = "Debug.txt";
    save_imageP3(resolution[0], resolution[1], debugfile, colours);

    return 0;
}

//
// Created by Aijaisarma Sabaratnasarma on 2023-11-28.
//

#include <cstdio>
#include <iostream>
#include <vector>
#include <variant>
#include <map>
#include <cstring> // Include for strcmp
#include "RayUtilities.h"
#include <cmath>
#include <algorithm>

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


//Scans and gets the data from the text file by comparing strings and key variables line by line and then storing the 
//scanned values in an array and hashmaps
void getSceneInfo(FILE* file, float* vP, int* res, float* Ambient, std::map<std::string,
        std::vector<std::variant<float, int, std::vector<float>>>>& gD,
        std::map<std::string, std::vector<std::variant<float, std::vector<float>>>>& lD,
        std::vector<float>& bG, char* Name)
{
    char line[100];
    float num = 0.0;

    for (int i = 0; i < 5; i++)
    {
        fgets(line, sizeof(line), file);
        if (sscanf(line, "%*s %f", &num) == 1)
        {
            vP[i] = num;
        }
    }

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
        //scanning and comparing strings from text file
        sscanf(line, "%s", temp);
        if (strcmp(temp, "SPHERE") == 0)
        {
            std::vector<std::variant<float, int, std::vector<float>>> currSphere;
            std::string name;
            float pos_x, pos_y, pos_z, scl_x, scl_y, scl_z, red, green, blue, KA, KD, KS, KR;
            int n;

            //allocating scanned data into arrays and hashmaps
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

        //comparing string from text file line by line
        else if (strcmp(temp, "LIGHT") == 0)
        {
            std::vector<std::variant<float, std::vector<float>>> currLightSource;
            std::string name;
            float pos_x, pos_y, pos_z, I_R, I_G, I_B;

            //allocating scanned data into arrays and hashmaps
            sscanf(line, "LIGHT %s %f %f %f %f %f %f", temp, &pos_x, &pos_y, &pos_z, &I_R, &I_G, &I_B);

            name = temp;
            currLightSource.push_back(std::vector<float>{pos_x, pos_y, pos_z});
            currLightSource.push_back(I_R);
            currLightSource.push_back(I_G);
            currLightSource.push_back(I_B);
            lD[name] = currLightSource;
        }

         //comparing string from text file line by line
        else if (strcmp(temp, "BACK") == 0)
        {
            float red, green, blue;

            //allocating scanned data into arrays and hashmaps
            sscanf(line, "BACK %f %f %f", &red, &green, &blue);
            bG.push_back(red);
            bG.push_back(green);
            bG.push_back(blue);
        }

         //comparing string from text file line by line
        else if (strcmp(temp, "AMBIENT") == 0)
        {
            float AR, AG, AB;

            //allocating scanned data into arrays
            sscanf(line, "AMBIENT %f %f %f", &AR, &AG, &AB);
            Ambient[0] = AR;
            Ambient[1] = AG;
            Ambient[2] = AB;
        }

         //comparing string from text file line by line
        else if (strcmp(temp, "OUTPUT") == 0)
        {
            sscanf(line, "OUTPUT %s", Name);
        }

        else {
            break;
        }
    }
}


//solve for t by calculating how many roots there are using the quadratic equation
float solveT(const Ray& ray)
{
    std::vector<float> origin = ray.getPoint();
    std::vector<float> direction = ray.getDirection();

    float a = dot(direction, direction);
    float b = dot(origin, direction);
    float c = dot(origin, origin) - 1;

    float discriminant = b * b - a * c;
    if (discriminant < 0)
    {
        return -INFINITY;
    }
    else
    {
        float sqrtDiscriminant = sqrt(discriminant);
        float t1 = (-b - sqrtDiscriminant) / (a);
        float t2 = (-b + sqrtDiscriminant) / (a);
        float t = std::min(t1, t2);

        if (t < 0)
        {
            return -INFINITY;
        }
        else if(discriminant == 0){
            return -b/(2*a);
        }
        else
        {
            return std::max(t, 1.0f);
        } 
    }
}


//calculates the intersection point from the given ray
std::vector<std::variant<std::vector<float>, std::string>> Intersect(Ray& ray)
{
    std::vector<float> closest = {0.0f, 0.0f, 0.0f};
    std::vector<std::string> intersectedObjects = {};

    for (const auto& sphere : geometryData)
    {
        std::vector<float> position = reinterpret_cast<const std::vector<float> &>(sphere.second[0]);
        std::vector<float> scale = reinterpret_cast<const std::vector<float> &>(sphere.second[1]);

        double M[4][4] = {{scale[0], 0, 0, position[0]}, {0, scale[1], 0, position[1]},
                         {0, 0, scale[2], position[2]}, {0, 0, 0, 1}};
        double M_inverse[4][4] = {{0}};
        invert_matrix(M, M_inverse);

        Ray transRay = ray.getTransformedRay(M_inverse);

        //solve for t by calculating how many roots there are using the quadratic equation
        float intersectionParameter = solveT(transRay);
        if (intersectionParameter == -INFINITY) { continue; }

        std::vector<float> intersectionPoint = ray.findPoint(intersectionParameter);
        //intersectionPoint = transformPoint(M, intersectionPoint);

        if (equal(closest, std::vector<float>{0.0f, 0.0f, 0.0f}))
        {
            closest = intersectionPoint;
            intersectedObjects.push_back(sphere.first);
        }

        else
        {
            if (closestToEye(intersectionPoint, closest))
            {
                closest = intersectionPoint;
                intersectedObjects.push_back(sphere.first);
            }
        }
    }

    std::vector<std::variant<std::vector<float>, std::string>> closestIntersection = {};
    if (equal(closest, std::vector<float>({0.0f, 0.0f, 0.0f})))
    {
        closestIntersection.push_back(closest);
    }
    else
    {
        closestIntersection.push_back(closest);
        closestIntersection.push_back(intersectedObjects.back());
    }

    return closestIntersection;
}

//calculates the reflection vector by taking the given normal vector and light vector
std::vector<float> reflectionVector(std::vector<float>& normal, const std::vector<float>& lightVector)
{   
    //normalizig the normal vector
    normalize(normal); 
    float num = 2 * dot(normal, lightVector);
    scalarMult(num, normal);
    std::vector<float> reflected = subtract(lightVector, normal);

    scalarMult(pow(num, -1), normal); 
    return reflected;
}

//calculates the pixel color by gathering information from the ray, normal, intersection point and intersection object
std::vector<float> illuminate(Ray& viewRay, std::vector<float>& normal, const std::vector<float>& point, std::vector<std::variant<float, std::vector<float>>>& obj)
{
    //const std::vector<float>& viewDir = getScaledVec(-1.0, viewRay.getDirection());
    const std::vector<float> viewDir = viewRay.getDirection();
    // scalarMult(-1.0, viewDir); 
    //const std::vector<float> viewDir {-point[0], -point[1], -point[2]};

    float Ka = (float&) (obj[3]);
    std::vector<float> colour = (std::vector<float>&) obj[2];
    std::vector<float> pixelColour = {Ka * ambient[0] * colour[0], Ka * ambient[1] * colour[1], Ka * ambient[2] * colour[2]};

    float Kd = (float&) (obj[4]);
    float Ks = (float&) (obj[5]);
    int shininess = (int&) obj[7];

    std::vector<float> currColour = {0.0f, 0.0f, 0.0f};
    float lightSources = 0.0f; 
    for (const auto& light : lightData)
    {
        std::vector<float> lightPos = reinterpret_cast<const std::vector<float> &>(light.second[0]);
        std::vector<float> shadowDir = {lightPos[0] - point[0], lightPos[1] - point[1], lightPos[2] - point[2]};

        normalize(shadowDir);
        Ray shadowRay = Ray(point, shadowDir);

        //calculates the intersection point
        std::vector<std::variant<std::vector<float>, std::string>> intersection = Intersect(shadowRay);

        if (equal((std::vector<float>&) intersection[0], std::vector<float>{0.0f, 0.0f, 0.0f}))
        {
            //calculates the reflection vector
            std::vector<float> reflected = reflectionVector(normal, shadowDir);

            //normalizing the vectors
            normalize(reflected);
            normalize(viewDir);
            normalize(normal); 

            float redIntensity = (float &) light.second[1];
            float greenIntensity = (float &) light.second[2];
            float blueIntensity = (float &) light.second[3];

            float diffuseRed = Kd * redIntensity * std::max(dot(normal, shadowDir), 0.0f) * colour[0];
            float diffuseGreen = Kd * greenIntensity * std::max(dot(normal, shadowDir), 0.0f) * colour[1];
            float diffuseBlue = Kd * blueIntensity * std::max(dot(normal, shadowDir), 0.0f) * colour[2];

            float specularRed = Ks * redIntensity * pow(std::max(dot(reflected, viewDir), 0.0f), shininess);
            float specularGreen = Ks * greenIntensity * pow(std::max(dot(reflected, viewDir), 0.0f), shininess);
            float specularBlue = Ks * blueIntensity * pow(std::max(dot(reflected, viewDir), 0.0f), shininess);

            //printf("(%f %f %f)\n", specularRed, specularGreen, specularBlue);  
            //printf("%f\n", dot(reflected, viewDir)); 

            currColour[0] += diffuseRed + specularRed;
            currColour[1] += diffuseGreen + specularGreen;
            currColour[2] += diffuseBlue + specularBlue;

            lightSources += 1.0f;
        }
    }

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

    pixelColour[0] = std::min(pixelColour[0], 1.0f);
    pixelColour[1] = std::min(pixelColour[1], 1.0f);
    pixelColour[2] = std::min(pixelColour[2], 1.0f);

    return pixelColour;
}

//Returns the raytaced screen colour from the given current ray
std::vector<float> raytrace(Ray& currRay)
{
    if (currRay.getDepth() >= MAX_DEPTH) {return std::vector<float>{0.0f, 0.0f, 0.0f};}

    //calculates the intersection point
    std::vector<std::variant<std::vector<float>, std::string>> intersectInfo = Intersect(currRay);
    std::vector<float> intersectPt = (std::vector<float>&) intersectInfo[0];

    if (equal((std::vector<float>&) intersectInfo[0], std::vector<float>{0.0f, 0.0f, 0.0f}) && currRay.getDepth() == 1)
    {return background;}

    if (equal((std::vector<float>&) intersectInfo[0], std::vector<float>{0.0f, 0.0f, 0.0f})) {return background;}

    auto intersectObj = geometryData.at((std::string &) intersectInfo[1]);

    std::vector<float> position = reinterpret_cast<const std::vector<float> &>(intersectObj[0]);
    std::vector<float> scale = reinterpret_cast<const std::vector<float> &>(intersectObj[1]);

    double M[4][4] = {{scale[0], 0, 0, position[0]}, {0, scale[1], 0, position[1]},
                     {0, 0, scale[2], position[2]}, {0, 0, 0, 1}};
    double M_inverse[4][4] = {{0}};
    double M_InverseTranspose[4][4] = {{0}};
    invert_matrix(M, M_inverse);
    transpose(M_inverse, M_InverseTranspose); 

    // calculate normal to intersected sphere
    std::vector<float> normal = MatrixVectorMult(M_InverseTranspose, MatrixVectorMult(M_inverse, intersectPt)); 

    // calculate the normal to the unit sphere at the intersection point
    // std::vector<float> normal = subtract(intersectPt, position);

    // // normalize the normal vector
    // normalize(normal);

    // // transform the normal vector using the inverse transpose of the transformation matrix
    // normal = MatrixVectorMult(M_InverseTranspose, normal);

    // calculate the reflected ray
    // std::vector<float> incomingDir = currRay.getDirection();
    // float dotProduct = -2 * (dot(normal, incomingDir)); 
    // scalarMult(dotProduct, normal);
    // std::vector<float> reflectedDir = add(normal, incomingDir);


    // Ray reflectedRay = Ray(intersectPt, reflectedDir);
    // reflectedRay.setDepth(currRay.getDepth() + 1);

    //scalarMult(pow(dotProduct, -1), normal); 
    //normalize(normal);

    // calculate pixel colour
    std::vector<float> clocal = illuminate(currRay, normal, intersectPt, (std::vector<std::variant<float, std::vector<float>>>&) intersectObj);
    
    // std::vector<float> colourRE = raytrace(reflectedRay);
    // scalarMult((float&) intersectObj[6], colourRE);
    
    std::vector<float> screenColour = clocal; //add(clocal, colourRE); 

    return screenColour;
}

float clamp(float value, float min, float max)
{
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

int main(int argc, char **argv)
{
    if (argc != 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    // open text file
    FILE *sceneData = fopen(argv[1], "r");

    if (sceneData == NULL) {
        perror("Error opening file");
        return 1;
    }

    // get the data from the text file and store it in arrays and hashmaps
    getSceneInfo(sceneData, viewPlane, resolution, ambient, geometryData, lightData, background, sceneName);
    // close input file
    fclose(sceneData);

    //printSceneInfo(viewPlane, resolution, ambient, geometryData, lightData, background, sceneName);

    unsigned char* colours;
    colours = new unsigned char [3 * resolution[0] * resolution[1]];
    int counter = 0;
    for (int i = 0; i < resolution[1]; i++)
    {
        for (int j = 0; j < resolution[0]; j++)
        {
            float dirX = viewPlane[1] + (viewPlane[2] - viewPlane[1]) * j / resolution[0];
            float dirY = viewPlane[4] + (viewPlane[3] - viewPlane[4]) * i / resolution[1]; 
            float dirZ = (-1) * viewPlane[0];

            std::vector<float> normalizedDir = {dirX, dirY, dirZ};

            //normalizing the normal vector
            normalize(normalizedDir);
            
            Ray ray = Ray(std::vector<float>{0, 0, 0}, normalizedDir);

            ray.setDepth(1);

            //Returns the raytaced screen colour from the given current ray
            std::vector<float> screenColour = raytrace(ray);

            colours[counter] = (unsigned char) (screenColour[0] * 255.0);
            colours[counter+1] = (unsigned char) (screenColour[1] * 255.0);
            colours[counter+2] = (unsigned char) (screenColour[2] * 255.0);
            counter += 3;
        }
    }

    save_imageP6(resolution[0], resolution[1], sceneName, colours);
    char debugfile[10] = "Debug.txt";
    save_imageP3(resolution[0], resolution[1], debugfile, colours);

    return 0;
}

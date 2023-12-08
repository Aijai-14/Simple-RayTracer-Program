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
        sscanf(line, "%s", temp);
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

        else if (strcmp(temp, "BACK") == 0)
        {
            float red, green, blue;
            sscanf(line, "BACK %f %f %f", &red, &green, &blue);
            bG.push_back(red);
            bG.push_back(green);
            bG.push_back(blue);
        }

        else if (strcmp(temp, "AMBIENT") == 0)
        {
            float AR, AG, AB;
            sscanf(line, "AMBIENT %f %f %f", &AR, &AG, &AB);
            Ambient[0] = AR;
            Ambient[1] = AG;
            Ambient[2] = AB;
        }

        else if (strcmp(temp, "OUTPUT") == 0)
        {
            sscanf(line, "OUTPUT %s", Name);
        }

        else {
            break;
        }
    }
}

float solveT(const Ray& ray)
{
    float a = dot(ray.getDirection(), ray.getDirection());
    float b = 2 * dot(ray.getPoint(), ray.getDirection());
    float c = dot(ray.getPoint(), ray.getPoint()) - 1;

    float discriminant = (float) pow(b, 2) - 4 * a * c;
    if (discriminant < 0)
    {
        return -INFINITY;
    }

    else if (discriminant == 0)
    {
        if ((-b / (2*a)) > 1) {return (-b / (2*a));}
        else {return -INFINITY;}
    }

    else
    {
        float t1 = (-b + sqrt(discriminant)) / (2*a);
        float t2 = (-b - sqrt(discriminant)) / (2*a);

        if (t1 < t2 && t1 > 1) {return t1;}
        else if (t2 < t1 && t2 > 1) {return t2;}
        else {return -INFINITY;}
    }
}

std::vector<std::variant<std::vector<float>, std::string>> Intersect(Ray& ray)
{
    std::vector<float> closest = {0.0f, 0.0f, 0.0f};
    std::vector<std::string> intersectedObjects = {};
    for (const auto& sphere : geometryData)
    {
        std::vector<float> position = reinterpret_cast<const std::vector<float> &>(sphere.second.at(0));
        std::vector<float> scale = reinterpret_cast<const std::vector<float> &>(sphere.second.at(1));

        double M[4][4] = {{scale[0], 0, 0, position[0]}, {0, scale[1], 0, position[1]},
                         {0, 0, scale[2], position[2]}, {0, 0, 0, 1}};
        double M_inverse[4][4] = {{0}};
        invert_matrix(M, M_inverse);

        Ray transRay = ray.getTransformedRay(M_inverse);

        float intersectionParameter = solveT(transRay);
        if (intersectionParameter == -INFINITY) { continue; }

        std::vector<float> intersectionPoint = ray.findPoint(intersectionParameter);

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

std::vector<float> reflectionVector(const std::vector<float>& normal, const std::vector<float>& lightVector)
{
    float num = 2 * dot(normal, lightVector);
    scalarMult(num, normal);
    std::vector<float> reflected = subtract(normal, lightVector);

    normalize(normal);
    return reflected;
}

std::vector<float> illuminate(Ray& viewRay, const std::vector<float>& normal, const std::vector<float>& point, std::vector<std::variant<float, std::vector<float>>>& obj)
{
    std::vector<float> viewDir = {((std::vector<float>) viewRay.getDirection())[0], ((std::vector<float>) viewRay.getDirection())[1],
                                  ((std::vector<float>) viewRay.getDirection())[2]};
    normalize(viewDir);
    normalize(normal);

    float Ka = (float&) (obj[3]);
    std::vector<float> colour = (std::vector<float>&) obj[2];
    std::vector<float> pixelColour = {Ka * ambient[0] * colour[0], Ka * ambient[1] * colour[1], Ka * ambient[2] * colour[2]};

    float Kd = (float&) (obj[4]);
    float Ks = (float&) (obj[5]);
    int shininess = (int&) obj[7];

    for (const auto& light : lightData)
    {
        std::vector<float> currColour = {0.0f, 0.0f, 0.0f};
        std::vector<float> lightPos = reinterpret_cast<const std::vector<float> &>(light.second[0]);
        std::vector<float> shadowDir = {lightPos[0] - point[0], lightPos[1] - point[1], lightPos[2] - point[2]};

        Ray shadowRay = Ray(point, shadowDir);
        std::vector<std::variant<std::vector<float>, std::string>> intersection = Intersect(shadowRay);

        if (equal((std::vector<float>&) intersection[0], std::vector<float>{0.0f, 0.0f, 0.0f}))
        {
            normalize(shadowDir);
            std::vector<float> reflected = reflectionVector(normal, shadowDir);

            float redIntensity = (float &) light.second[1];
            float greenIntensity = (float &) light.second[2];
            float blueIntensity = (float &) light.second[3];

            float diffuseRed = Kd * redIntensity * (dot(normal, shadowDir)) * colour[0];
            float specularRed = (float) (Ks * redIntensity * pow(dot(reflected, viewDir), shininess));
            float diffuseGreen = Kd * greenIntensity * (dot(normal, shadowDir)) * colour[1];
            float specularGreen = (float) (Ks * greenIntensity * pow(dot(reflected, viewDir), shininess));
            float diffuseBlue = Kd * blueIntensity * (dot(normal, shadowDir)) * colour[2];
            float specularBlue = (float) (Ks * blueIntensity * pow(dot(reflected, viewDir), shininess));

            currColour[0] = diffuseRed + specularRed;
            currColour[1] = diffuseGreen + specularGreen;
            currColour[2] = diffuseBlue + specularBlue;
        }

        for (int i = 0; i < 3; i++)
        {
            pixelColour[i] += currColour[i];
        }
    }

    if (pixelColour[0] > 1.0) {pixelColour[0] = 1.0;}
    else if (pixelColour[1] > 1.0) {pixelColour[1] = 1.0;}
    else if (pixelColour[2] > 1.0) {pixelColour[2] = 1.0;}

    return pixelColour;
}

std::vector<float> raytrace(Ray& currRay)
{
    if (currRay.getDepth() >= MAX_DEPTH) {return std::vector<float>{0.0f, 0.0f, 0.0f};}

    std::vector<std::variant<std::vector<float>, std::string>> intersectInfo = Intersect(currRay);
    std::vector<float> intersectPt = (std::vector<float>&) intersectInfo[0];

    if (equal((std::vector<float>&) intersectInfo[0], std::vector<float>{0.0f, 0.0f, 0.0f}) || currRay.getDepth() == 1)
    {return background;}

    std::vector<std::variant<float, std::vector<float>>> intersectObj =
            reinterpret_cast<const std::vector<std::variant<float, std::vector<float>>> &>(geometryData.at((std::string &) intersectInfo[1]));

    // calculate normal to intersected sphere
    std::vector<float> center = (std::vector<float>&) intersectObj[0];
    std::vector<float> normal = {2 * (intersectPt[0] - center[0]), 2 * (intersectPt[1] - center[1]), 2 * (intersectPt[2] - center[2])};

    // calculate the reflected ray
    std::vector<float> incomingDir = currRay.getDirection();
    scalarMult(-2 * (dot(normal, incomingDir)), normal);
    std::vector<float> reflectedDir = add(normal, incomingDir);
    Ray reflectedRay = Ray(intersectPt, reflectedDir);
    reflectedRay.setDepth(currRay.getDepth() + 1);

    // calculate the refracted ray

    // calculate pixel colour
    std::vector<float> clocal = illuminate(currRay, normal, intersectPt, intersectObj);
    std::vector<float> colourRE = raytrace(reflectedRay);

    scalarMult((float&) intersectObj[6], colourRE);
    std::vector<float> screenColour = add(clocal, colourRE);
    return screenColour;
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
    getSceneInfo(sceneData, viewPlane, resolution, ambient,geometryData, lightData, background, sceneName);
    // close input file
    fclose(sceneData);

    unsigned char* colours;
    colours = new unsigned char [3 * resolution[0] * resolution[1]];
    int counter = 0;
    for (int i = 1; i < resolution[0] + 1; i++)
    {
        for (int j = 1; j < resolution[1] + 1; j++)
        {
            float dirX = viewPlane[2] * (((float) (2 * i) / (float) resolution[0]) - 1);
            float dirY = viewPlane[4] * (((float) (2 * j) / (float) resolution[1]) - 1);
            float dirZ = (-1) * viewPlane[0];
            Ray ray = Ray(std::vector<float>{0.000001, 0.000001, 0.000001}, std::vector<float>{dirX, dirY, dirZ});

            ray.setDepth(1);

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

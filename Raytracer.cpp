//
// Created by Aijaisarma Sabaratnasarma on 2023-11-28.
//

#include <cstdio>
#include <iostream>
#include <vector>
#include <variant>
#include <map>
#include <cstring> // Include for strcmp
#include "RayClass.h"
#include <cmath>

const int MAX_DEPTH = 4;
std::map<std::string, std::vector<std::variant<float, int, std::vector<float>>>> geometryData;
std::map<std::string, std::vector<std::variant<float, std::vector<float>>>> lightData;
std::vector<float> background;

void getSceneInfo(FILE* file, float* vP, int* res, float* Ambient, std::map<std::string,
        std::vector<std::variant<float, int, std::vector<float>>>>& gD,
        std::map<std::string, std::vector<std::variant<float, std::vector<float>>>>& lD,
        std::vector<float>& bG, std::string& Name)
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
            char outputFile[20];
            sscanf(line, "OUTPUT %s", outputFile);
            Name = outputFile;
        }

        else {
            break;
        }
    }
}

void printSceneData(float* vP, int* res, std::map<std::string,
        std::vector<std::variant<float, int, std::vector<float>>>>& gD,
        std::map<std::string, std::vector<std::variant<float, std::vector<float>>>>& lD,
        std::vector<float>& bG, float* Ambient, std::string& Name)
{
    std::cout << Name << std::endl;
    std::cout << "Plane Parameters: ";
    for (int i = 0; i < 5; i++)
    {
        std::cout << vP[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Resolution: ";
    for (int i = 0; i < 2; i++)
    {
        std::cout << res[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Ambient Intensities: ";
    for (int i = 0; i < 3; i++)
    {
        std::cout << Ambient[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Background Colour: ";
    for (float colour : bG)
    {
        std::cout << colour << " ";
    }
    std::cout << std::endl;

    for (const auto& entry : gD) {
        const std::string& name = entry.first;
        const std::vector<std::variant<float, int, std::vector<float>>>& values = entry.second;

        // Print key
        std::cout << "Sphere: " << name << std::endl;

        // Loop through the values in the vector
        for (const auto& value : values) {
            if (std::holds_alternative<float>(value)) {
                std::cout << "  Light Coefficients (adsr): " << std::get<float>(value) << std::endl;
            }
            else if (std::holds_alternative<int>(value)) {
                std::cout << "  Shininess: " << std::get<int>(value) << std::endl;
            }
            else if (std::holds_alternative<std::vector<float>>(value)) {
                const auto& innerVector = std::get<std::vector<float>>(value);
                std::cout << "  Position/Scale/Colour: ";
                for (float innerValue : innerVector) {
                    std::cout << innerValue << " ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;  // Add a separator between objects
    }

    for (const auto& entry : lD) {
        const std::string& name = entry.first;
        const std::vector<std::variant<float, std::vector<float>>>& values = entry.second;

        // Print key
        std::cout << "Light: " << name << std::endl;

        // Loop through the values in the vector
        for (const auto& value : values) {
            if (std::holds_alternative<float>(value)) {
                std::cout << "  Light Intensities (rgb): " << std::get<float>(value) << std::endl;
            }

            else if (std::holds_alternative<std::vector<float>>(value)) {
                const auto& innerVector = std::get<std::vector<float>>(value);
                std::cout << "  Position: ";
                for (float innerValue : innerVector) {
                    std::cout << innerValue << " ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;  // Add a separator between objects
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
        return (-b / (2*a));
    }

    else
    {
        float t1 = (-b + sqrt(discriminant)) / (2*a);
        float t2 = (-b - sqrt(discriminant)) / (2*a);

        if (t1 < t2) {return t1;}
        else {return t2;}
    }
}

std::vector<float> Intersect(Ray& ray, std::map<std::string, std::vector<std::variant<float, int, std::vector<float>>>>& spheres)
{
    std::vector<float> closest = {0.0f, 0.0f, 0.0f};
    for (const auto& sphere : spheres)
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
        }

        else
        {
            if (closestToEye(intersectionPoint, closest))
            {
                closest = intersectionPoint;
            }
        }
    }

    return closest;
}


//std::vector<float> raytrace(Ray& currRay)
//{
//    if (currRay.getDepth() > MAX_DEPTH) {return std::vector<float>{0.0f, 0.0f, 0.0f};}
//    std::vector<float> intersectPt = Intersect(currRay, geometryData);
//    if (equal(intersectPt, std::vector<float>{0.0f, 0.0f, 0.0f}) && currRay.getDepth() == 1) {return background;}
//
//
//
//
//
//
//}

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

    // create array for near plane parameters, array for resolution, hashmaps for sphere data and light data, vector for
    // background colour, array for ambient light intensities, and string for output file name.
    float viewPlane [5] = {};
    int resolution [2] = {};
    float ambient [3] = {};
//    std::map<std::string, std::vector<std::variant<float, int, std::vector<float>>>> geometryData;
//    std::map<std::string, std::vector<std::variant<float, std::vector<float>>>> lightData;
//    std::vector<float> background;
    std::string sceneName;

    // get the data from the text file and store it in arrays and hashmaps
    getSceneInfo(sceneData, viewPlane, resolution, ambient,geometryData, lightData, background, sceneName);

    // close input file
    fclose(sceneData);

    // print scene data to test extraction
    //printSceneData(viewPlane, resolution, geometryData, lightData, background, ambient, sceneName);

//    Intersect(geometryData);

//    std::vector<float> colours(3 * resolution[0] * resolution[1], 0.0f);
//    for (int i = 0; i < resolution[0]; i++)
//    {
//        for (int j = 0; j < resolution[1]; j++)
//        {
//            float dirX = viewPlane[2] * (((float) (2 * i) / (float) resolution[0]) - 1);
//            float dirY = viewPlane[4] * (((float) (2 * j) / (float) resolution[1]) - 1);
//            float dirZ = (-1) * viewPlane[0];
//            Ray ray = Ray(std::vector<float>{0, 0, 0}, std::vector<float>{dirX, dirY, dirZ});
//
//            ray.setDepth(1);
//
////            tracedRay = raytrace(ray);
//
//        }
//    }

    return 0;
}

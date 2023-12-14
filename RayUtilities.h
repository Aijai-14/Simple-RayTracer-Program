#include <utility>
#include <cmath>

//
// Created by Aijaisarma Sabaratnasarma on 2023-11-29.
//

#ifndef SIMPLE_RAYTRACER_PROGRAM_RAYCLASS_H
#define SIMPLE_RAYTRACER_PROGRAM_RAYCLASS_H

// Output in P6 format, a binary file containing:
// P6
// ncolumns nrows
// Max colour value
// colours in binary format thus unreadable
void save_imageP6(int Width, int Height, char* fname,unsigned char* pixels) {
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname,Width,Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3,Width,fp);
    }

    fclose(fp);
}

void printSceneInfo(float* vP, int* res, float* Ambient, 
                    std::map<std::string, std::vector<std::variant<float, int, std::vector<float>>>>& gD,
                    std::map<std::string, std::vector<std::variant<float, std::vector<float>>>>& lD,
                    std::vector<float>& bG, char* Name)
{
    std::cout << "Viewport: ";
    for (int i = 0; i < 5; i++) {
        std::cout << vP[i] << " ";
    }
    std::cout << "\nResolution: " << res[0] << " " << res[1] << "\n";
    std::cout << "Ambient: " << Ambient[0] << " " << Ambient[1] << " " << Ambient[2] << "\n";
    std::cout << "Background: " << bG[0] << " " << bG[1] << " " << bG[2] << "\n";
    std::cout << "Output: " << Name << "\n";

    std::cout << "Spheres:\n";
    for (const auto& sphere : gD) {
        std::cout << "Name: " << sphere.first << "\n";
        for (const auto& data : sphere.second) {
            if (std::holds_alternative<int>(data)) {
                std::cout << "n: " << std::get<int>(data) << "\n";
            } else if (std::holds_alternative<float>(data)) {
                std::cout << "KA, KD, KS, KR: " << std::get<float>(data) << "\n";
            } else if (std::holds_alternative<std::vector<float>>(data)) {
                std::vector<float> vec = std::get<std::vector<float>>(data);
                std::cout << "Vector: ";
                for (float f : vec) {
                    std::cout << f << " ";
                }
                std::cout << "\n";
            }
        }
    }

    std::cout << "Lights:\n";
    for (const auto& light : lD) {
        std::cout << "Name: " << light.first << "\n";
        for (const auto& data : light.second) {
            if (std::holds_alternative<float>(data)) {
                std::cout << "I_R, I_G, I_B: " << std::get<float>(data) << "\n";
            } else if (std::holds_alternative<std::vector<float>>(data)) {
                std::vector<float> vec = std::get<std::vector<float>>(data);
                std::cout << "Position: ";
                for (float f : vec) {
                    std::cout << f << " ";
                }
                std::cout << "\n";
            }
        }
    }
}

// Output in P3 format, a text file containing:
// P3
// ncolumns nrows
// Max colour value (for us, and usually 255)
// r1 g1 b1 r2 g2 b2 .....
void save_imageP3(int Width, int Height, char* fname,unsigned char* pixels) {
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname,Width,Height);
    fp = fopen(fname,"w");
    if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
    }
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    int k = 0 ;
    for(int j = 0; j < Height; j++) {

        for( int i = 0 ; i < Width; i++)
        {
            fprintf(fp," %d %d %d", pixels[k],pixels[k+1],pixels[k+2]) ;
            k = k + 3 ;
        }
        fprintf(fp,"\n") ;
    }
    fclose(fp);
}

std::vector<float> add(std::vector<float> vec1, std::vector<float> vec2)
{
    return std::vector<float>{vec1[0] + vec2[0], vec1[1] + vec2[1], vec1[2] + vec2[2]};
}

std::vector<float> subtract(std::vector<float> vec1, std::vector<float> vec2)
{
    return std::vector<float>{vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2]};
}

void scalarMult(double constant, std::vector<float> vec)
{
    for (int i = 0; i < 3; i++)
    {
        vec[i] = constant * vec[i];
    }
}

std::vector<float> getScaledVec(double constant, std::vector<float> vec)
{
    std::vector<float> scaledVec = vec;
    for (int i = 0; i < 3; i++)
    {
        scaledVec[i] = constant * vec[i];
    }
    return scaledVec;
}

void normalize(std::vector<float> vec)
{
    auto magnitude = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
    auto inverse_magnitude = pow(magnitude, -1);
    scalarMult(inverse_magnitude, vec);
}

bool closestToEye(std::vector<float> vec1, std::vector<float> vec2) {

    auto dis1 = (float) sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
    auto dis2 = (float) sqrt(pow(vec2[0], 2) + pow(vec2[1], 2) + pow(vec2[2], 2));

    if (dis1 < dis2) {return true;}
    else {return false;}
}

bool equal(std::vector<float> vec1, std::vector<float> vec2)
{
    for (int i = 0; i < 3; i++)
    {
        if (vec1[i] != vec2[i])
        {
            return false;
        }
    }

    return true;
}

float dot(const double *vec1, std::vector<float> vec2)
{

    float product = 0.0f;
    for (int i = 0; i < 4; i++)
    {
        product += (float) vec1[i] * vec2[i];
    }
    return product;
}

float dot(std::vector<float> vec1, std::vector<float> vec2)
{

    float product = 0.0f;
    for (int i = 0; i < 3; i++)
    {
        product += vec1[i] * vec2[i];
    }
    return product;
}

std::vector<float> MatrixVectorMult(double M[][4], const std::vector<float>& vec)
{
    std::vector<float> output = {0.0f, 0.0f, 0.0f, 0.0f};

    if (vec.size() == 3) 
    {
        std::vector<float> homogenousVec = vec; 
        homogenousVec.push_back(1.0f);  

        for (int i = 0; i < 4; i++)
        {
            output[i] = dot(M[i], homogenousVec);
        }
    }

    else 
    {
        for (int i = 0; i < 4; i++)
        {
            output[i] = dot(M[i], vec);
        }
    }

    return output;
}

class Ray {
    private:
        std::vector<float> point;
        std::vector<float> direction;
        int depth = 0;

    public:
        Ray(std::vector<float> pt, std::vector<float> dir)
            {
                    point = std::move(pt);
                    direction = std::move(dir);
            }

    public:
        void setDepth(int d)
        {
            depth = d;
        }

    public:
        int getDepth() const
        {
            return depth;
        }

public:
        Ray getTransformedRay(double M[][4])
        {
            std::vector<float> homogenousPos = point;
            homogenousPos.push_back(1.0f);
            std::vector<float> homogenousDir = direction;
            homogenousDir.push_back(0.0f);

            std::vector<float> homogenousTransPos = MatrixVectorMult(M, homogenousPos);
            std::vector<float> homogenousTransDir = MatrixVectorMult(M, homogenousDir);

            std::vector<float> TransPos = {homogenousTransPos[0], homogenousTransPos[1], homogenousTransPos[2]};
            std::vector<float> TransDir = {homogenousTransDir[0], homogenousTransDir[1], homogenousTransDir[2]};
            //normalize(TransDir);

            return {TransPos, TransDir};
        }

public:
    const std::vector<float> &getPoint() const {
        return point;
    }

public:
    const std::vector<float> &getDirection() const {
        return direction;
    }

public:
    std::vector<float> findPoint(float t)
    {
       std::vector<float> pt = {0.0f, 0.0f, 0.0f};
       for (int i = 0; i < 3; i++)
       {
           pt[i] = point[i] + (t * direction[i]);
       }

        return pt;
    }
};

/*------------------------------------------------------*/
/*
 * Invert a 4x4 matrix.  Changed slightly from Richard Carling's code
 * in "Graphics Gems I".
 */

#define SMALL_NUMBER    1.e-8

void transpose(double M[4][4], double M_transpose[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            M_transpose[j][i] = M[i][j];
        }
    }
}

//double det2x2( a, b, c, d)
//   double a, b, c, d;
double det2x2( double a, double b, double c, double d)
{
    double ans;
    ans = a * d - b * c;
    return ans;
}

//double det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3 )
//   double a1, a2, a3, b1, b2, b3, c1, c2, c3;
double det3x3( double a1, double a2, double a3, double b1, double b2, double b3, double c1,
               double c2, double c3 )
{
    double ans;

    ans = a1 * det2x2( b2, b3, c2, c3 )
          - b1 * det2x2( a2, a3, c2, c3 )
          + c1 * det2x2( a2, a3, b2, b3 );
    return ans;
}

double det4x4( double m[4][4] )
{
    double ans;
    double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, d4;


    a1 = m[0][0]; b1 = m[0][1];
    c1 = m[0][2]; d1 = m[0][3];

    a2 = m[1][0]; b2 = m[1][1];
    c2 = m[1][2]; d2 = m[1][3];

    a3 = m[2][0]; b3 = m[2][1];
    c3 = m[2][2]; d3 = m[2][3];

    a4 = m[3][0]; b4 = m[3][1];
    c4 = m[3][2]; d4 = m[3][3];

    ans = a1 * det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4)
          - b1 * det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4)
          + c1 * det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4)
          - d1 * det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);
    return ans;
}

void adjoint( double in[4][4], double out[4][4] ) {
    double a1, a2, a3, a4, b1, b2, b3, b4;
    double c1, c2, c3, c4, d1, d2, d3, d4;


    a1 = in[0][0]; b1 = in[0][1];
    c1 = in[0][2]; d1 = in[0][3];

    a2 = in[1][0]; b2 = in[1][1];
    c2 = in[1][2]; d2 = in[1][3];

    a3 = in[2][0]; b3 = in[2][1];
    c3 = in[2][2]; d3 = in[2][3];

    a4 = in[3][0]; b4 = in[3][1];
    c4 = in[3][2]; d4 = in[3][3];

    out[0][0]  =   det3x3( b2, b3, b4, c2, c3, c4, d2, d3, d4);
    out[1][0]  = - det3x3( a2, a3, a4, c2, c3, c4, d2, d3, d4);
    out[2][0]  =   det3x3( a2, a3, a4, b2, b3, b4, d2, d3, d4);
    out[3][0]  = - det3x3( a2, a3, a4, b2, b3, b4, c2, c3, c4);

    out[0][1]  = - det3x3( b1, b3, b4, c1, c3, c4, d1, d3, d4);
    out[1][1]  =   det3x3( a1, a3, a4, c1, c3, c4, d1, d3, d4);
    out[2][1]  = - det3x3( a1, a3, a4, b1, b3, b4, d1, d3, d4);
    out[3][1]  =   det3x3( a1, a3, a4, b1, b3, b4, c1, c3, c4);

    out[0][2]  =   det3x3( b1, b2, b4, c1, c2, c4, d1, d2, d4);
    out[1][2]  = - det3x3( a1, a2, a4, c1, c2, c4, d1, d2, d4);
    out[2][2]  =   det3x3( a1, a2, a4, b1, b2, b4, d1, d2, d4);
    out[3][2]  = - det3x3( a1, a2, a4, b1, b2, b4, c1, c2, c4);

    out[0][3]  = - det3x3( b1, b2, b3, c1, c2, c3, d1, d2, d3);
    out[1][3]  =   det3x3( a1, a2, a3, c1, c2, c3, d1, d2, d3);
    out[2][3]  = - det3x3( a1, a2, a3, b1, b2, b3, d1, d2, d3);
    out[3][3]  =   det3x3( a1, a2, a3, b1, b2, b3, c1, c2, c3);
}

void invert_matrix (double A[4][4], double Ainv[4][4]) {
    int i, j;
    double det ;

    adjoint( A, Ainv );

    det = det4x4( A );

    if ( fabs( det ) < SMALL_NUMBER) {
        fprintf(stderr,"invert_matrix: matrix is singular!");
        return;
    }

    for (i=0; i<4; i++)
        for(j=0; j<4; j++)
            Ainv[i][j] = Ainv[i][j] / det;
}


#endif //SIMPLE_RAYTRACER_PROGRAM_RAYCLASS_H

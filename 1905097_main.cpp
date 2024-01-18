#include<bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;


// Point structure
struct Point{
    double x, y, z, n;

    // Constructor
    Point(double x, double y, double z){
        this->x = x; 
        this->y = y; 
        this->z = z;
        this->n = 1;
    }

    Point(){
        this->x = 0; 
        this->y = 0; 
        this->z = 0;
        this->n = 1;
    }

    // Copy constructor
    Point(const Point &p){
        this->x = p.x; 
        this->y = p.y; 
        this->z = p.z;
        this->n = p.n;
    }
};


// Add two points
Point addTwoPoints(Point a, Point b){
    return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}


// Subtract two points
Point subtractTwoPoints(Point a, Point b){
    return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}


// Multiply point with scalar
Point multiplyPointWithScalar(Point a, double s){
    return Point(a.x * s, a.y * s, a.z * s);
}


// Divide point with scalar
Point dividePointWithScalar(Point a, double s){
    return Point(a.x / s, a.y / s, a.z / s);
}


// Dot product of two points
double dotProduct(Point a, Point b){
    double res = a.x * b.x + a.y * b.y + a.z * b.z;
    return res;
}


// Cross product of two points
Point crossProduct(Point a, Point b){
    Point c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}


// Normalize a point
Point normalizePoint(Point a){
    double length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return Point(a.x / length, a.y / length, a.z / length);
}


// Random number generator
static unsigned long int g_seed = 1;
inline int random(){
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}


// Rotate a point
Point rotatePoint(Point p, Point axis, double angle){
    angle = angle * acos(-1) / 180.0;
    axis = normalizePoint(axis);

    Point p1 = addTwoPoints(multiplyPointWithScalar(p, cos(angle)), multiplyPointWithScalar(crossProduct(axis, p), sin(angle)));
    p1 = addTwoPoints(p1, multiplyPointWithScalar(axis, dotProduct(axis, p) * (1 - cos(angle))));
    return p1;
}


// Triangle structure
struct Triangle{
    Point a, b, c;
    int colorR, colorG, colorB;

    // Constructor
    Triangle(Point a, Point b, Point c){
        this->a = a;
        this->b = b;
        this->c = c;
    }

    Triangle(){
        this->a = Point();
        this->b = Point();
        this->c = Point();
    }

    // Assign Color
    void setColor(){
        colorR = random() % 256;
        colorG = random() % 256;
        colorB = random() % 256;
    }

    // Sort points based on y
    void sortBasedOnY(){
        if(a.y > b.y)
            swap(a, b);
        if(a.y > c.y)
            swap(a, c);
        if(b.y > c.y)
            swap(b, c);
    }
};


// Matrix structure
struct Matrix{
    double matrix[4][4];


    // Constructor
    Matrix(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                matrix[i][j] = 0;
        }
    }


    // Copy constructor
    Matrix(const Matrix &m){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                matrix[i][j] = m.matrix[i][j];
        }
    }


    // Set identity matrix
    void setIdentityMatrix(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                if(i == j)
                    matrix[i][j] = 1;
                else
                    matrix[i][j] = 0;
            }    
        }
    }


    // Translate matrix
    void translate(Point p){
        matrix[0][3] = p.x;
        matrix[1][3] = p.y;
        matrix[2][3] = p.z;
    }


    // Scale matrix
    void scale(Point p){
        matrix[0][0] = p.x;
        matrix[1][1] = p.y;
        matrix[2][2] = p.z;
    }


    // Rotate matrix using Rodrigues' rotation formula
    void rotate(double angle, Point axis){
        axis = normalizePoint(axis);
        setIdentityMatrix();

        Point i(1.0, 0.0, 0.0), j(0.0, 1.0, 0.0), k(0.0, 0.0, 1.0);
        i = rotatePoint(i, axis, angle);
        j = rotatePoint(j, axis, angle);
        k = rotatePoint(k, axis, angle);

        matrix[0][0] = i.x;
        matrix[0][1] = j.x;
        matrix[0][2] = k.x;
        matrix[1][0] = i.y;
        matrix[1][1] = j.y;
        matrix[1][2] = k.y;
        matrix[2][0] = i.z;
        matrix[2][1] = j.z;
        matrix[2][2] = k.z;
    }
};


// Multiply two matrices
Matrix multiplyTwoMatrix(Matrix a, Matrix b){
    Matrix m;

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            m.matrix[i][j] = 0;
            for(int k = 0; k < 4; k++)
                m.matrix[i][j] += a.matrix[i][k] * b.matrix[k][j];
        }
    }
    return m;
}


// Multiply matrix with point
Point multiplyMatrixWithPoint(Matrix m, Point p){
    Point p1;
    p1.x = m.matrix[0][0] * p.x + m.matrix[0][1] * p.y + m.matrix[0][2] * p.z + m.matrix[0][3] * p.n;
    p1.y = m.matrix[1][0] * p.x + m.matrix[1][1] * p.y + m.matrix[1][2] * p.z + m.matrix[1][3] * p.n;
    p1.z = m.matrix[2][0] * p.x + m.matrix[2][1] * p.y + m.matrix[2][2] * p.z + m.matrix[2][3] * p.n;
    p1.n = m.matrix[3][0] * p.x + m.matrix[3][1] * p.y + m.matrix[3][2] * p.z + m.matrix[3][3] * p.n;
    return Point(p1.x / p1.n, p1.y / p1.n, p1.z / p1.n);
}


// Multiply matrix with triangle
Triangle multiplyMatrixWithTriangle(Matrix m, Triangle t){
    Triangle t1;
    t1.a = multiplyMatrixWithPoint(m, t.a);
    t1.b = multiplyMatrixWithPoint(m, t.b);
    t1.c = multiplyMatrixWithPoint(m, t.c);
    return t1;
}


// Main function
int main()
{

    // ********* Stage 1 ***********
    // Input and output files
    ifstream in("scene.txt");
    ofstream out("stage1.txt");

    // Input variables
    Point eye, look, up;
    double fovY, aspectRatio, near, far;

    // Read input from file
    in >> eye.x >> eye.y >> eye.z;
    in >> look.x >> look.y >> look.z;
    in >> up.x >> up.y >> up.z;
    in >> fovY >> aspectRatio >> near >> far;

    // Initialize stack and matrix
    Matrix mat;
    mat.setIdentityMatrix();
    stack<Matrix> st;
    st.push(mat);
    string command;
    Point a, b, c;

    // Read commands from file
    while(in >> command){

        // Read triangle and transform using matrix and output to file
        if(command == "triangle"){
            in >> a.x >> a.y >> a.z;
            in >> b.x >> b.y >> b.z;
            in >> c.x >> c.y >> c.z;
            Triangle triangle(a, b, c);
            triangle = multiplyMatrixWithTriangle(st.top(), triangle);
            out << setprecision(7) << fixed << triangle.a.x << " " << triangle.a.y << " " << triangle.a.z << " \n" << triangle.b.x << " " << triangle.b.y << " " << triangle.b.z << " \n" << triangle.c.x << " " << triangle.c.y << " " << triangle.c.z << " \n" << endl;
        }

        // Translate the top matrix of stack
        else if(command == "translate"){
            Point p;
            in >> p.x >> p.y >> p.z;
            mat.setIdentityMatrix();
            mat.translate(p);
            st.top() = multiplyTwoMatrix(st.top(), mat);
        }

        // Scale the top matrix of stack
        else if(command == "scale"){
            Point p;
            in >> p.x >> p.y >> p.z;
            mat.setIdentityMatrix();
            mat.scale(p);
            st.top() = multiplyTwoMatrix(st.top(), mat);
        }

        // Rotate the top matrix of stack
        else if(command == "rotate"){
            double angle;
            Point axis;
            in >> angle >> axis.x >> axis.y >> axis.z;
            mat.setIdentityMatrix();
            mat.rotate(angle, axis);
            st.top() = multiplyTwoMatrix(st.top(), mat);
        }

        // Push the top matrix of stack
        else if(command == "push"){
            st.push(st.top());
        }

        // Pop the top matrix of stack
        else if(command == "pop"){
            st.pop();
        }

        // End the program
        else if(command == "end"){
            break;
        }

    }

    // Close input and output files
    in.close();
    out.close();


    // ********* Stage 2 ***********
    // Input and output files
    in.open("stage1.txt");
    out.open("stage2.txt");

    // Mutually perpendicular unit vectors
    Point l, r, u;

    // Calculate l, r, u
    l = subtractTwoPoints(look, eye);
    l = normalizePoint(l);

    r = crossProduct(l, up);
    r = normalizePoint(r);

    u = crossProduct(r, l);
    u = normalizePoint(u);

    // Move the eye to origin
    Matrix transformMatrix;
    transformMatrix.setIdentityMatrix();
    Point eyeNeg = multiplyPointWithScalar(eye, -1);
    transformMatrix.translate(eyeNeg);

    // Calculate the View Transformation matrix
    Matrix R;
    R.setIdentityMatrix();
    R.matrix[0][0] = r.x;
    R.matrix[0][1] = r.y;
    R.matrix[0][2] = r.z;
    R.matrix[1][0] = u.x;
    R.matrix[1][1] = u.y;
    R.matrix[1][2] = u.z;
    R.matrix[2][0] = -l.x;
    R.matrix[2][1] = -l.y;
    R.matrix[2][2] = -l.z;

    Matrix v = multiplyTwoMatrix(R, transformMatrix);

    // Read triangles from file, transform using matrix and output to file
    while(in >> a.x){
        in >> a.y >> a.z >> b.x >> b.y >> b.z >> c.x >> c.y >> c.z;
        Triangle triangle(a, b, c);
        triangle = multiplyMatrixWithTriangle(v, triangle);
        out << setprecision(7) << fixed << triangle.a.x << " " << triangle.a.y << " " << triangle.a.z << " \n" << triangle.b.x << " " << triangle.b.y << " " << triangle.b.z << " \n" << triangle.c.x << " " << triangle.c.y << " " << triangle.c.z << " \n" << endl;
    }

    // Close input and output files
    in.close();
    out.close();


    // ********* Stage 3 ***********
    // Input and output files
    in.open("stage2.txt");
    out.open("stage3.txt");

    // Calculate the Projection Transformation matrix
    double fovX = fovY * aspectRatio;
    fovX = fovX * acos(-1) / 180.0;
    fovY = fovY * acos(-1) / 180.0;
    double t1 = near * tan(fovY / 2);
    double r2 = near * tan(fovX / 2);

    Matrix projectionMatrix;
    projectionMatrix.matrix[0][0] = near / r2;
    projectionMatrix.matrix[1][1] = near / t1;
    projectionMatrix.matrix[2][2] = -(far + near) / (far - near);
    projectionMatrix.matrix[2][3] = -(2.0 * far * near) / (far - near);
    projectionMatrix.matrix[3][2] = -1.0;

    // Read triangles from file, transform using Projection Transformation Matrix and output to file
    while(in >> a.x){
        in >> a.y >> a.z >> b.x >> b.y >> b.z >> c.x >> c.y >> c.z;
        Triangle triangle(a, b, c);
        triangle = multiplyMatrixWithTriangle(projectionMatrix, triangle);
        out << setprecision(7) << fixed << triangle.a.x << " " << triangle.a.y << " " << triangle.a.z << " \n" << triangle.b.x << " " << triangle.b.y << " " << triangle.b.z << " \n" << triangle.c.x << " " << triangle.c.y << " " << triangle.c.z << " \n" << endl;
    }

    // Close input and output files
    in.close();
    out.close();


    // ********* Stage 4 ***********
    // Input and output files
    in.open("stage3.txt");
    out.open("z_buffer.txt");
    ifstream in2("config.txt");
    int screenWidth, screenHeight;

    // Read screen width and height from file
    in2 >> screenWidth >> screenHeight;

    // Initialize Screen Parameters
    double topY = 1.0, bottomY = -1.0, leftX = -1.0, rightX = 1.0;
    double dx = (rightX - leftX) / screenWidth;
    double dy = (topY - bottomY) / screenHeight;
    double topPixelY = topY - dy / 2.0;
    double leftPixelX = leftX + dx / 2.0;
    double zMin = -1.0, zMax = 1.0;

    // Initialize z-buffer
    vector<vector<double>> zBuffer(screenHeight, vector<double>(screenWidth));
    for(int i = 0; i < screenHeight; i++){
        for(int j = 0; j < screenWidth; j++)
            zBuffer[i][j] = zMax;
    }

    // Initialize image
    bitmap_image image(screenWidth, screenHeight);
    image.set_all_channels(0, 0, 0);

    // Read triangles from file, Update the z buffer values
    while(in >> a.x){
        in >> a.y >> a.z >> b.x >> b.y >> b.z >> c.x >> c.y >> c.z;
        Triangle triangle(a, b, c);

        // Sort the points of triangle based on y
        triangle.sortBasedOnY();

        // Set color of triangle
        triangle.setColor();

        // Calculate the yMin and yMax
        double yMin = triangle.a.y;
        double yMax = triangle.c.y;
        yMin = max(yMin, bottomY);
        yMax = min(yMax, topY);

        // If three points are parallel to y axis, continue
        if(yMin == yMax){
           continue; 
        }

        // Calculate the top and bottom scanline
        int topScanline = (int)ceil((topPixelY - yMax) / dy);
        int bottomScanline = (int)floor((topPixelY - yMin) / dy);

        // Calculate the points for each row
        for(int rowNo = topScanline; rowNo <= bottomScanline; rowNo++){

            // Calculate the y value of the row
            double yVal = topPixelY - rowNo * dy;
            Point p1, p2;

            // Calculate the intersection points of the row with the triangle
            if(yVal >= triangle.a.y && yVal <= triangle.b.y && triangle.a.y != triangle.b.y){
                p1 = Point(triangle.a.x + (triangle.b.x - triangle.a.x) * (yVal - triangle.a.y) / (triangle.b.y - triangle.a.y), yVal, triangle.a.z + (triangle.b.z - triangle.a.z) * (yVal - triangle.a.y) / (triangle.b.y - triangle.a.y));
            }
            else if(yVal >= triangle.b.y && yVal <= triangle.c.y && triangle.b.y != triangle.c.y){
                p1 = Point(triangle.b.x + (triangle.c.x - triangle.b.x) * (yVal - triangle.b.y) / (triangle.c.y - triangle.b.y), yVal, triangle.b.z + (triangle.c.z - triangle.b.z) * (yVal - triangle.b.y) / (triangle.c.y - triangle.b.y));
            }

            p2 = Point(triangle.a.x + (triangle.c.x - triangle.a.x) * (yVal - triangle.a.y) / (triangle.c.y - triangle.a.y), yVal, triangle.a.z + (triangle.c.z - triangle.a.z) * (yVal - triangle.a.y) / (triangle.c.y - triangle.a.y));

            // Sort the points based on x
            if(p1.x > p2.x)
                swap(p1, p2);

            // Calculate the xMin and xMax
            double xMin = max(p1.x, leftX);
            double xMax = min(p2.x, rightX);

            // Calculate the left and right column  
            int leftColumn = (int)round((xMin - leftPixelX) / dx);
            int rightColumn = (int)round((xMax - leftPixelX) / dx);

            // Calculate the z value of each pixel and update the z buffer
            for(int colNo = leftColumn; colNo <= rightColumn; colNo++){
                double xVal = leftPixelX + colNo * dx;
                double zVal = p2.z;
                if(p1.x != p2.x)
                    zVal += (p1.z - p2.z) * (p2.x - xVal) / (p2.x - p1.x);

                if(zVal < zBuffer[rowNo][colNo] && zVal > zMin){
                    zBuffer[rowNo][colNo] = zVal;
                    image.set_pixel(colNo, rowNo, triangle.colorR, triangle.colorG, triangle.colorB);
                }
            }

        }       
    }

    // Save the image as bmp file
    image.save_image("out.bmp");

    // Output the z buffer to file
    for(int i = 0; i < screenHeight; i++){
        for(int j = 0; j < screenWidth; j++){
            if(zBuffer[i][j] < zMax){
                out << setprecision(6) << fixed << zBuffer[i][j] << "\t";
            }    
        }    
        out << endl;
    }

    // Close input and output files
    in.close();
    in2.close();
    out.close();

    // Clear the z buffer
    zBuffer.clear();
    zBuffer.resize(0);

    return 0;
}
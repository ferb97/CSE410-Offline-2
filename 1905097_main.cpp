#include<bits/stdc++.h>
using namespace std;

struct Point{
    double x, y, z, n;

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

    Point(const Point &p){
        this->x = p.x; 
        this->y = p.y; 
        this->z = p.z;
        this->n = p.n;
    }
};

Point addTwoPoints(Point a, Point b){
    return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

Point subtractTwoPoints(Point a, Point b){
    return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}

Point multiplyPointWithScalar(Point a, double s){
    return Point(a.x * s, a.y * s, a.z * s);
}

Point dividePointWithScalar(Point a, double s){
    return Point(a.x / s, a.y / s, a.z / s);
}

double dotProduct(Point a, Point b){
    double res = a.x * b.x + a.y * b.y + a.z * b.z;
    return res;
}

Point crossProduct(Point a, Point b){
    Point c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

Point normalizePoint(Point a){
    double length = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    return Point(a.x / length, a.y / length, a.z / length);
}

struct Triangle{
    Point a, b, c;

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
};

struct Matrix{
    double matrix[4][4];

    Matrix(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                matrix[i][j] = 0;
        }
    }

    Matrix(const Matrix &m){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                matrix[i][j] = m.matrix[i][j];
        }
    }

    void print(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                cout << matrix[i][j] << " ";
            cout << endl;
        }
    }

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

    void translate(Point p){
        matrix[0][3] = p.x;
        matrix[1][3] = p.y;
        matrix[2][3] = p.z;
    }

    void scale(Point p){
        matrix[0][0] = p.x;
        matrix[1][1] = p.y;
        matrix[2][2] = p.z;
    }

    void rotate(double angle, Point axis){
        angle = angle * acos(-1) / 180.0;
        double c = cos(angle);
        double s = sin(angle);
        double t = 1 - c;

        axis = normalizePoint(axis);
        setIdentityMatrix();

        matrix[0][0] = t * axis.x * axis.x + c;
        matrix[0][1] = t * axis.x * axis.y - s * axis.z;
        matrix[0][2] = t * axis.x * axis.z + s * axis.y;
        matrix[1][0] = t * axis.x * axis.y + s * axis.z;
        matrix[1][1] = t * axis.y * axis.y + c;
        matrix[1][2] = t * axis.y * axis.z - s * axis.x;
        matrix[2][0] = t * axis.x * axis.z - s * axis.y;
        matrix[2][1] = t * axis.y * axis.z + s * axis.x;
        matrix[2][2] = t * axis.z * axis.z + c;
    }
};

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

Point multiplyMatrixWithPoint(Matrix m, Point p){
    Point p1;
    p1.x = m.matrix[0][0] * p.x + m.matrix[0][1] * p.y + m.matrix[0][2] * p.z + m.matrix[0][3] * p.n;
    p1.y = m.matrix[1][0] * p.x + m.matrix[1][1] * p.y + m.matrix[1][2] * p.z + m.matrix[1][3] * p.n;
    p1.z = m.matrix[2][0] * p.x + m.matrix[2][1] * p.y + m.matrix[2][2] * p.z + m.matrix[2][3] * p.n;
    p1.n = m.matrix[3][0] * p.x + m.matrix[3][1] * p.y + m.matrix[3][2] * p.z + m.matrix[3][3] * p.n;
    return Point(p1.x / p1.n, p1.y / p1.n, p1.z / p1.n);
}

Triangle multiplyMatrixWithTriangle(Matrix m, Triangle t){
    Triangle t1;
    t1.a = multiplyMatrixWithPoint(m, t.a);
    t1.b = multiplyMatrixWithPoint(m, t.b);
    t1.c = multiplyMatrixWithPoint(m, t.c);
    return t1;
}

Matrix getIdentityMatrix(){
    Matrix m;
    for(int i = 0; i < 4; i++)
        m.matrix[i][i] = 1;
    return m;
}

int main()
{
    ifstream in("scene.txt");
    ofstream out("stage1.txt");

    Point eye, look, up;
    double fovY, aspectRatio, near, far;

    in >> eye.x >> eye.y >> eye.z;
    in >> look.x >> look.y >> look.z;
    in >> up.x >> up.y >> up.z;
    in >> fovY >> aspectRatio >> near >> far;

    Matrix mat;
    mat.setIdentityMatrix();
    stack<Matrix> st;
    st.push(mat);
    string command;
    Point a, b, c;

    while(in >> command){

        if(command == "triangle"){
            in >> a.x >> a.y >> a.z;
            in >> b.x >> b.y >> b.z;
            in >> c.x >> c.y >> c.z;
            Triangle triangle(a, b, c);
            triangle = multiplyMatrixWithTriangle(st.top(), triangle);
            out << setprecision(7) << fixed << triangle.a.x << " " << triangle.a.y << " " << triangle.a.z << "\n" << triangle.b.x << " " << triangle.b.y << " " << triangle.b.z << "\n" << triangle.c.x << " " << triangle.c.y << " " << triangle.c.z << endl << endl;
        }

        else if(command == "translate"){
            Point p;
            in >> p.x >> p.y >> p.z;
            mat.setIdentityMatrix();
            mat.translate(p);
            st.top() = multiplyTwoMatrix(st.top(), mat);
        }

        else if(command == "scale"){
            Point p;
            in >> p.x >> p.y >> p.z;
            mat.setIdentityMatrix();
            mat.scale(p);
            st.top() = multiplyTwoMatrix(st.top(), mat);
        }

        else if(command == "rotate"){
            double angle;
            Point axis;
            in >> angle >> axis.x >> axis.y >> axis.z;
            mat.setIdentityMatrix();
            mat.rotate(angle, axis);
            st.top() = multiplyTwoMatrix(st.top(), mat);
        }

        else if(command == "push"){
            st.push(st.top());
        }

        else if(command == "pop"){
            st.pop();
        }

        else if(command == "end"){
            break;
        }

    }

    in.close();
    out.close();

    in.open("stage1.txt");
    out.open("stage2.txt");
    Point l, r, u;

    l = subtractTwoPoints(look, eye);
    l = normalizePoint(l);

    r = crossProduct(l, up);
    r = normalizePoint(r);

    u = crossProduct(r, l);
    u = normalizePoint(u);

    Matrix transformMatrix;
    transformMatrix.setIdentityMatrix();
    Point eyeNeg = multiplyPointWithScalar(eye, -1);
    transformMatrix.translate(eyeNeg);

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

    while(in >> a.x){
        in >> a.y >> a.z >> b.x >> b.y >> b.z >> c.x >> c.y >> c.z;
        Triangle triangle(a, b, c);
        triangle = multiplyMatrixWithTriangle(v, triangle);
        out << setprecision(7) << fixed << triangle.a.x << " " << triangle.a.y << " " << triangle.a.z << "\n" << triangle.b.x << " " << triangle.b.y << " " << triangle.b.z << "\n" << triangle.c.x << " " << triangle.c.y << " " << triangle.c.z << endl << endl;
    }

    in.close();
    out.close();

    in.open("stage2.txt");
    out.open("stage3.txt");

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

    while(in >> a.x){
        in >> a.y >> a.z >> b.x >> b.y >> b.z >> c.x >> c.y >> c.z;
        Triangle triangle(a, b, c);
        triangle = multiplyMatrixWithTriangle(projectionMatrix, triangle);
        out << setprecision(7) << fixed << triangle.a.x << " " << triangle.a.y << " " << triangle.a.z << "\n" << triangle.b.x << " " << triangle.b.y << " " << triangle.b.z << "\n" << triangle.c.x << " " << triangle.c.y << " " << triangle.c.z << endl << endl;
    }

    in.close();
    out.close();

    return 0;
}
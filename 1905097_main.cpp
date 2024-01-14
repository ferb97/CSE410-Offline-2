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
    Point c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}

Point subtractTwoPoints(Point a, Point b){
    Point c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}

Point multiplyPointWithScalar(Point a, double s){
    Point c;
    c.x = a.x * s;
    c.y = a.y * s;
    c.z = a.z * s;
    return c;
}

Point dividePointWithScalar(Point a, double s){
    Point c;
    c.x = a.x / s;
    c.y = a.y / s;
    c.z = a.z / s;
    return c;
}

double dotProduct(Point a, Point b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
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
    Point c;
    c.x = a.x / length;
    c.y = a.y / length;
    c.z = a.z / length;
    return c;
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
    double mat[4][4];

    Matrix(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                mat[i][j] = 0;
        }
    }

    Matrix(const Matrix &m){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                mat[i][j] = m.mat[i][j];
        }
    }

    void print(){
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++)
                cout << mat[i][j] << " ";
            cout << endl;
        }
    }

    void setIdentityMatrix(){
        for(int i = 0; i < 4; i++)
            mat[i][i] = 1;
    }

    void translate(Point p){
        mat[0][3] = p.x;
        mat[1][3] = p.y;
        mat[2][3] = p.z;
    }

    void scale(Point p){
        mat[0][0] = p.x;
        mat[1][1] = p.y;
        mat[2][2] = p.z;
    }

    void rotate(double angle, Point axis){
        angle = angle * acos(-1) / 180.0;
        double c = cos(angle);
        double s = sin(angle);
        double t = 1 - c;
        axis = normalizePoint(axis);
        setIdentityMatrix();
        mat[0][0] = t * axis.x * axis.x + c;
        mat[0][1] = t * axis.x * axis.y - s * axis.z;
        mat[0][2] = t * axis.x * axis.z + s * axis.y;
        mat[1][0] = t * axis.x * axis.y + s * axis.z;
        mat[1][1] = t * axis.y * axis.y + c;
        mat[1][2] = t * axis.y * axis.z - s * axis.x;
        mat[2][0] = t * axis.x * axis.z - s * axis.y;
        mat[2][1] = t * axis.y * axis.z + s * axis.x;
        mat[2][2] = t * axis.z * axis.z + c;
    }
};

Matrix multiplyTwoMatrix(Matrix a, Matrix b){
    Matrix c;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            c.mat[i][j] = 0;
            for(int k = 0; k < 4; k++)
                c.mat[i][j] += a.mat[i][k] * b.mat[k][j];
        }
    }
    return c;
}

Point multiplyMatrixWithPoint(Matrix m, Point p){
    Point c;
    c.x = m.mat[0][0] * p.x + m.mat[0][1] * p.y + m.mat[0][2] * p.z + m.mat[0][3] * p.n;
    c.y = m.mat[1][0] * p.x + m.mat[1][1] * p.y + m.mat[1][2] * p.z + m.mat[1][3] * p.n;
    c.z = m.mat[2][0] * p.x + m.mat[2][1] * p.y + m.mat[2][2] * p.z + m.mat[2][3] * p.n;
    c.n = m.mat[3][0] * p.x + m.mat[3][1] * p.y + m.mat[3][2] * p.z + m.mat[3][3] * p.n;
    return c;
}

Matrix getIdentityMatrix(){
    Matrix c;
    for(int i = 0; i < 4; i++)
        c.mat[i][i] = 1;
    return c;
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

    Matrix m;
    m.setIdentityMatrix();
    stack<Matrix> s;
    s.push(m);
    string command;

    while(in >> command){
        if(command == "triangle"){
            Point a, b, c;
            in >> a.x >> a.y >> a.z;
            in >> b.x >> b.y >> b.z;
            in >> c.x >> c.y >> c.z;
            a.n = b.n = c.n = 1;
            a = multiplyMatrixWithPoint(s.top(), a);
            b = multiplyMatrixWithPoint(s.top(), b);
            c = multiplyMatrixWithPoint(s.top(), c);
            out << setprecision(7) << fixed << a.x << " " << a.y << " " << a.z << "\n" << b.x << " " << b.y << " " << b.z << "\n" << c.x << " " << c.y << " " << c.z << endl << endl;
        }
        else if(command == "translate"){
            Point p;
            in >> p.x >> p.y >> p.z;
            Matrix m;
            m.setIdentityMatrix();
            m.translate(p);
            s.top() = multiplyTwoMatrix(s.top(), m);
        }
        else if(command == "scale"){
            Point p;
            in >> p.x >> p.y >> p.z;
            Matrix m;
            m.setIdentityMatrix();
            m.scale(p);
            s.top() = multiplyTwoMatrix(s.top(), m);
        }
        else if(command == "rotate"){
            double angle;
            Point axis;
            in >> angle >> axis.x >> axis.y >> axis.z;
            Matrix m;
            m.setIdentityMatrix();
            m.rotate(angle, axis);
            s.top() = multiplyTwoMatrix(s.top(), m);
        }
        else if(command == "push"){
            s.push(s.top());
        }
        else if(command == "pop"){
            s.pop();
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

    Matrix t;
    t.setIdentityMatrix();
    Point eyeNeg = multiplyPointWithScalar(eye, -1);
    t.translate(eyeNeg);

    Matrix r1;
    r1.setIdentityMatrix();
    r1.mat[0][0] = r.x;
    r1.mat[0][1] = r.y;
    r1.mat[0][2] = r.z;
    r1.mat[1][0] = u.x;
    r1.mat[1][1] = u.y;
    r1.mat[1][2] = u.z;
    r1.mat[2][0] = -l.x;
    r1.mat[2][1] = -l.y;
    r1.mat[2][2] = -l.z;

    Matrix v = multiplyTwoMatrix(r1, t);
    Point a1, b1, c1;

    while(in >> a1.x){
        in >> a1.y >> a1.z >> b1.x >> b1.y >> b1.z >> c1.x >> c1.y >> c1.z;
        a1.n = b1.n = c1.n = 1;
        a1 = multiplyMatrixWithPoint(v, a1);
        b1 = multiplyMatrixWithPoint(v, b1);
        c1 = multiplyMatrixWithPoint(v, c1);
        out << setprecision(7) << fixed << a1.x << " " << a1.y << " " << a1.z << "\n" << b1.x << " " << b1.y << " " << b1.z << "\n" << c1.x << " " << c1.y << " " << c1.z << endl << endl;
    }

    in.close();
    out.close();
    
    return 0;
}
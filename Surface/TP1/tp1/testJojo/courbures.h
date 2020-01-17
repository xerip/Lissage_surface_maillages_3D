#ifndef COURBURES_H
#define COURBURES_H
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/ArrayKernel.hh>
#include <OpenMesh/Core/Geometry/Vector11T.hh>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include "mainwindow.h"

class MyQuad
{
    double _coefs[5] ; // a_0 x^2 + a1 xy + a2 y^2 + a3 x + a4 y + a5
    Eigen::AngleAxisd _r ;
    Eigen::Translation3d _t ;
public:
    MyQuad(double *data,
           const Eigen::AngleAxisd &r = Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)),
           const Eigen::Translation3d &t = Eigen::Translation3d(0,0,0))
        : _r(r), _t(t)
    {
        for (int i=0; i<5; i++)
            _coefs[i] = data[i] ;
    }
    MyQuad(const Eigen::VectorXd &v,
           const Eigen::AngleAxisd &r = Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)),
           const Eigen::Translation3d &t = Eigen::Translation3d(0,0,0))
        : _r(r), _t(t)
    {
        for (int i=0; i<5; i++)
            _coefs[i] = v[i] ;
    }
    MyQuad(double a0=0, double a1=0, double a2=0, double a3=0, double a4=0,
           const Eigen::AngleAxisd &r = Eigen::AngleAxisd(0, Eigen::Vector3d(0,0,1)),
           const Eigen::Translation3d &t = Eigen::Translation3d(0,0,0))
        : _r(r), _t(t)
    {
        _coefs[0] = a0 ;
        _coefs[1] = a1 ;
        _coefs[2] = a2 ;
        _coefs[3] = a3 ;
        _coefs[4] = a4 ;
    }
    MyQuad(const MyQuad &q) : _r(q._r), _t(q._t)
    {
        for (int i=0; i<5; i++)
            _coefs[i] = q._coefs[i] ;
    }

    double & operator[] (int i) {return _coefs[i] ; }

    double quad_fun (double x, double y)
    {
        return _coefs[0] * x*x + _coefs[1] * x*y + _coefs[2] * y*y + _coefs[3] * x + _coefs[4] * y  ;
    }

    double quad_fun (const OpenMesh::Vec2d &v)
    {

        return quad_fun(v[0], v[1]) ;
    }

};

template <typename T>
class MyStats
{
private:
    std::vector<T> _distrib ;
public:
    MyStats () {} ;

    void push_back (T data)
    {
        _distrib.push_back(data) ;
    }

    T min ()
    {
        T tmp = _distrib.at(0) ;
        for (int i=1 ; i<_distrib.size(); i++)
        {
            if (_distrib.at(i) < tmp)
                tmp = _distrib.at(i) ;
        }
        return tmp ;
    }

    T max ()
    {
        T tmp = _distrib.at(0) ;
        for (int i=1 ; i<_distrib.size(); i++)
        {
            if (_distrib.at(i) > tmp)
                tmp = _distrib.at(i) ;
        }
        return tmp ;
    }

    T mean ()
    {
        T acc(0) ;
        std::cout << "acc : " << acc << std::endl;
        if (_distrib.size() > 0)
        {
            for(int i=0; i<_distrib.size(); i++)
            {
                acc += _distrib.at(i) ;
            }
            return acc/(_distrib.size()) ;
        }
        else
            return acc ;
    }

    T stdev ()
    {
        T m = mean() ;
        T acc(0), tmp ;
        if (_distrib.size() > 0)
        {
            for(int i=0; i<size(_distrib); i++)
            {
                tmp = _distrib.at(i) - m ;
                acc += tmp * tmp ;
            }
            return sqrt(acc/(_distrib.size())) ;
        }
        else
            return acc ;
    }

    T stdev (T m)
    {
        T acc, tmp ;
        if (_distrib.size() > 0)
        {
            for(int i=0; i<_distrib.size(); i++)
            {
                tmp = _distrib.at(i) - m ;
                acc += tmp * tmp ;
            }
            return sqrt(acc/(_distrib.size())) ;
        }
        else
            return acc ;
    }
};

class Courbures
{
private:
    MyMesh &_mesh ;

    void trianguler(MyMesh::VertexHandle sommets, int size);

public:
    OpenMesh::VPropHandleT<double>    vprop_K;
    OpenMesh::VPropHandleT<double>    vprop_H;

    Courbures(MyMesh &mesh) : _mesh(mesh) {}

    void set_fixed_colors() ;
    void normales_locales() ;
    std::vector<MyMesh::VertexHandle> get_two_neighborhood(MyMesh::VertexHandle vh);
    MyQuad fit_quad(MyMesh::VertexHandle vh) ;
    int draw_quad(MyMesh::VertexHandle vh, MyMesh *mesh);
    void compute_KH() ;
    void set_K_colors() ;
};

#endif // COURBURES_H

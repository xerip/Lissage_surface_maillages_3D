#include <iostream>
#include "courbures.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>

void Courbures::set_fixed_colors()
{
    MyMesh::VertexIter v_it ;
    MyMesh::Point  p ;
    MyMesh::VertexHandle vh ;

    int i = 0 ;
    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        i++ ;
        if (i % 2)
            _mesh.set_color(*v_it, MyMesh::Color(0, 255, 0)) ;
        else
            _mesh.set_color(*v_it, MyMesh::Color(0, 0, 255)) ;
    }
}

void Courbures::normales_locales()
{
    _mesh.request_vertex_normals() ;

    int i ;
    MyMesh::Normal n ;

    for (MyMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        n = MyMesh::Normal(0,0,0) ;
        i = 0 ;
        for(MyMesh::VertexFaceIter vf_it = _mesh.vf_iter(*v_it);vf_it.is_valid();++vf_it)
        {
            i++ ;
            n += _mesh.calc_face_normal(*vf_it) ;
        }
        if (i!=0)
        {
            n /= i ;
            n = n / norm(n) ;
        }
        _mesh.set_normal(*v_it, n) ;
    }
}

std::vector<MyMesh::VertexHandle> Courbures::get_two_neighborhood(const MyMesh::VertexHandle vh)
{
    OpenMesh::VPropHandleT<bool>    vprop_flag;
    _mesh.add_property(vprop_flag,  "vprop_flag");

    // Initialisation
    for (MyMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        _mesh.property(vprop_flag, *v_it) = false ;
    }
    // Circulateur sur le premier cercle
    std::vector<MyMesh::VertexHandle> neigh, neigh2 ;

	_mesh.property(vprop_flag, vh) = true ;
    for(MyMesh::VertexVertexIter vv_it = _mesh.vv_iter(vh);vv_it.is_valid();++vv_it)
    {
        neigh.push_back(*vv_it) ; // ajout du point à la liste
        _mesh.property(vprop_flag, *vv_it) = true ;
    }

    // Parcours du premier cercle et ajout du second cercle par circulateurs
    for (int i=0; i<neigh.size(); i++)
    {
        MyMesh::VertexHandle vh = neigh.at(i) ;
        for(MyMesh::VertexVertexIter vv_it = _mesh.vv_iter(vh);vv_it.is_valid();++vv_it)
        {
            if (!_mesh.property(vprop_flag, *vv_it)) // sommet non encore rencontré
            neigh2.push_back(*vv_it) ; // ajout du point à la liste
            _mesh.property(vprop_flag, *vv_it) = true ;
        }
    }

    // Concaténation des deux cercles
    neigh.insert(neigh.end(), neigh2.begin(), neigh2.end()) ;

    _mesh.remove_property(vprop_flag);

    return neigh ;
}

MyQuad Courbures::fit_quad(MyMesh::VertexHandle vh)
{
    MyQuad q ;

    std::vector<MyMesh::VertexHandle> neigh = get_two_neighborhood(vh) ;

//    std::cout << "-- Neighborhood" << std::endl ;
    for (int i = 0; i< neigh.size(); i++)
    {
        MyMesh::Point p = _mesh.point(neigh.at(i)) ;
//        std::cout << p[0] << ", " << p[1] << ", " << p[2] << " ; " << std::endl ;
    }

    // Calcul de la matrice de changement de base
    MyMesh::Normal n = _mesh.normal(vh) ;
    MyMesh::Point p = _mesh.point(vh);
//    std::cout << "-- normale" << std::endl ;
//    std::cout << n[0] << ", " << n[1] << ", " << n[2] << std::endl ;
//    std::cout << "-- point" << std::endl ;
//    std::cout << p[0] << ", " << p[1] << ", " << p[2] << std::endl ;

    Eigen::Vector3d ne(n[0], n[1], n[2]), Oz(0,0,1), axis;
    Eigen::Vector3d p_e(p[0], p[1], p[2]), pi_e ;

    axis = ne.cross(Oz) ;
    double sina = axis.norm(), cosa = ne.dot(Oz), angle ;
    if (sina >= 0)
        angle = acos(cosa) ;
    else
        angle = -acos(cosa) ;
    axis = axis.normalized() ;

    Eigen::AngleAxisd r(angle, axis) ;
//    std::cout << "-- rotation" << std:: endl ;
//    std::cout << r.matrix()(0,0) << ", " << r.matrix()(0, 1) << ", " << r.matrix()(0,2) << "; " ;
//    std::cout << r.matrix()(1,0) << ", " << r.matrix()(1, 1) << ", " << r.matrix()(1,2) << "; " ;
//    std::cout << r.matrix()(2,0) << ", " << r.matrix()(2, 1) << ", " << r.matrix()(2,2) << std::endl ;
    Eigen::Translation3d t(-p_e[0], -p_e[1], -p_e[2]) ;
    Eigen::Transform<double, 3, Eigen::Affine> ch_base = r * t ;


    // Calcul de la matrice / vecteur de moindres carrés linéaires (Eigen)

    if (neigh.size() >= 5)
    {
        int n(neigh.size()) ;
        Eigen::MatrixXd A(n,5);
        Eigen::VectorXd B(n);

//        std::cout << "-- Après changement de base" << std::endl ;
        for(int i=0; i<neigh.size(); i++)
        {
            MyMesh::Point p = _mesh.point(neigh.at(i)) ;
            pi_e << p[0], p[1], p[2] ;
            pi_e = ch_base * pi_e ; // Application du changement de base
//            std::cout << pi_e[0] << ", " << pi_e[1] << ", " << pi_e[2] << std::endl ;

            // TODO : initialisation de A / B
            A(i, 0) = pi_e(0) * pi_e(0);    // x carré
            A(i, 1) = pi_e(1) * pi_e(1);    // y carré
            A(i, 2) = pi_e(0) * pi_e(1);    // x * y
            A(i, 3) = pi_e(0);              // x
            A(i, 4) = pi_e(1);              // y

            B(i) = - pi_e(2);
        }

        // Résolution aux moindres carrés par SVD
        Eigen::VectorXd coef(5) ;
        coef = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B) ;
        MyQuad q(coef[0], coef[1], coef[2], coef[3], coef[4], r, t) ;
//        std::cout << "-- quad : " << std::endl ;
//        std::cout << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << ", " << q[4] << ", " << q[5] << ", " << std::endl ;
        return q ;
    }
    else
    {
        std::cout << "Quad fitting : not enough neighbors" ;
        throw "Quad fitting : not enough neighbors" ;
    }
}

void Courbures::draw_quad(MyMesh::VertexHandle vh, MyMesh *mesh)
{
    MyQuad q = fit_quad(vh);

    std::vector<MyMesh::VertexHandle> neigh = get_two_neighborhood(vh) ;

    // Calcul de la matrice de changement de base
    MyMesh::Normal n = _mesh.normal(vh) ;
    MyMesh::Point p = _mesh.point(vh);

    Eigen::Vector3d ne(n[0], n[1], n[2]), Oz(0,0,1), axis;
    Eigen::Vector3d p_e(p[0], p[1], p[2]), pi_e ;

    axis = ne.cross(Oz) ;
    double sina = axis.norm(), cosa = ne.dot(Oz), angle ;
    if (sina >= 0)
        angle = acos(cosa) ;
    else
        angle = -acos(cosa) ;
    axis = axis.normalized() ;

    Eigen::AngleAxisd r(angle, axis) ;;
    Eigen::Translation3d t(-p_e[0], -p_e[1], -p_e[2]) ;
    Eigen::Transform<double, 3, Eigen::Affine> ch_base = r * t ;


    // Calcul de la matrice / vecteur de moindres carrés linéaires (Eigen)

    if (neigh.size() >= 5)
    {
        MyMesh::VertexHandle sommets[neigh.size()];
        for(int i=0; i<neigh.size(); i++)
        {
            MyMesh::VertexHandle vh = neigh.at(i);

            mesh->data(vh).thickness = 5;
            mesh->set_color(vh, MyMesh::Color(0, 255, 0));

            MyMesh::Point p1 = _mesh.point(neigh.at(i)) ;
            pi_e << p1[0], p1[1], p1[2] ;
            pi_e = ch_base * pi_e ; // Application du changement de base

            Eigen::Vector3d tmpp(pi_e[0], pi_e[1], q.quad_fun(pi_e[0],pi_e[1]));
//            tmpp = ch_base * tmpp;

            //ajout sommet quadrique
            sommets[i] = mesh->add_vertex(MyMesh::Point(tmpp[0], tmpp[1], tmpp[2]));

        }

        std::vector<MyMesh::VertexHandle> uneNouvelleFace;

        for (int i  = 0; i < neigh.size(); i++){
            for (int j  = 0; j < neigh.size(); j++){
                uneNouvelleFace.clear();
                uneNouvelleFace.push_back(sommets[i]);
                uneNouvelleFace.push_back(sommets[(j)]);
                uneNouvelleFace.push_back(sommets[(j+1)%(neigh.size()-1)]);
                mesh->add_face(uneNouvelleFace);

                MyMesh::FaceHandle fh = mesh->face_handle(mesh->n_faces()-1);
                mesh->request_face_colors();
                mesh->set_color(fh,MyMesh::Color(255,100,0));
            }
        }
    }
    else
    {
        std::cout << "Quad fitting : not enough neighbors" ;
        throw "Quad fitting : not enough neighbors" ;
    }
}


void Courbures::compute_KH()
{
    MyMesh::VertexIter v_it ;
    MyQuad  q ;
    MyMesh::VertexHandle vh ;

    _mesh.add_property(vprop_K,  "vprop_K");
    _mesh.add_property(vprop_H,  "vprop_H");

    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        vh = *v_it ;
        q = fit_quad(vh) ;
        _mesh.property(vprop_K, vh) = 4 * q[0] * q[2] - q[1]*q[1] ; // K = det (K_P) = 4 a_0 a_2 - a_1 ^ 2
        _mesh.property(vprop_H, vh) = q[0] + q[2] ; // K = Tr (K_P) = a_0 + a_2
    }
}

void Courbures::set_K_colors()
{
    MyMesh::VertexIter v_it ;
    MyStats<double> dist ;
    OpenMesh::Vec3uc tmp_color ;
    double r = 255, g = 255, b = 255 ;

    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        dist.push_back(_mesh.property(vprop_K, *v_it));
    }
    double m = dist.mean(), sigma = dist.stdev(m) ;
    double tmp, lambda ;
    std::cout << "K min : " << dist.min() << " | max : " << dist.max() << std::endl ;
    std::cout << "K mean : " << m << " | sigma : " << sigma << std::endl;

    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        tmp  = _mesh.property(vprop_K, *v_it) - m ;
		// Carte de couleur de tmp (pour les valeurs comprises entre
		// m-sigma / m+sigma

        // TODO

	   if (tmp < 0)
        {
            b = 255;
            r = 255 * (-sigma+tmp)/-sigma;
            g = r;
        }
        else if (tmp > 0)
        {
            r = 255;
            b = 255 * (sigma-tmp)/sigma;
            g = b;
        }
        else {
            r = 255; g = 255; b = 255 ;
        }

		tmp_color = OpenMesh::Vec3uc(r,g,b) ;
        _mesh.set_color(*v_it, MyMesh::Color(tmp_color)) ;
    }
}

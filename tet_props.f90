module tet_props
  implicit none

  private


  public :: tetrahedron_dihedral_angles, tetrahedron_edge_length
  public :: tetrahedron_volume

contains

  function r8_acos ( c )

    !*****************************************************************************80
    !
    !! R8_ACOS computes the arc cosine function, with argument truncation.
    !
    !  Discussion:
    !
    !    If you call your system ACOS routine with an input argument that is
    !    even slightly outside the range [-1.0, 1.0 ], you may get an unpleasant
    !    surprise (I did).
    !
    !    This routine simply truncates arguments outside the range.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 October 2012
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) C, the argument.
    !
    !    Output, real ( kind = 8 ) R8_ACOS, an angle whose cosine is C.
    !
    implicit none

    real ( kind = 8 ) c
    real ( kind = 8 ) c2
    real ( kind = 8 ) r8_acos

    c2 = c
    c2 = max ( c2, -1.0D+00 )
    c2 = min ( c2, +1.0D+00 )

    r8_acos = acos ( c2 )

    return
  end function r8_acos

  subroutine r8_swap ( x, y )

    !*****************************************************************************80
    !
    !! R8_SWAP swaps two R8's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 December 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
    !    Y have been interchanged.
    !
    implicit none

    real ( kind = 8 ) x
    real ( kind = 8 ) y
    real ( kind = 8 ) z

    z = x
    x = y
    y = z

    return
  end subroutine r8_swap

  function r8mat_det_4d ( a )

    !*****************************************************************************80
    !
    !! R8MAT_DET_4D computes the determinant of a 4 by 4 matrix.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    16 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
    !
    !    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
    !
    implicit none

    real ( kind = 8 ) a(4,4)
    real ( kind = 8 ) r8mat_det_4d

    r8mat_det_4d = &
         a(1,1) * ( &
         a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
         - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
         + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
         - a(1,2) * ( &
         a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
         - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
         + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
         + a(1,3) * ( &
         a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
         - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
         + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
         - a(1,4) * ( &
         a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
         - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
         + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

    return
  end function r8mat_det_4d

  subroutine r8mat_solve ( n, rhs_num, a, info )

    !*****************************************************************************80
    !
    !! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    29 August 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the order of the matrix.
    !
    !    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.
    !    RHS_NUM must be at least 0.
    !
    !    Input/output, real ( kind = 8 ) A(N,N+rhs_num), contains in rows and
    !    columns 1 to N the coefficient matrix, and in columns N+1 through
    !    N+rhs_num, the right hand sides.  On output, the coefficient matrix
    !    area has been destroyed, while the right hand sides have
    !    been overwritten with the corresponding solutions.
    !
    !    Output, integer ( kind = 4 ) INFO, singularity flag.
    !    0, the matrix was not singular, the solutions were computed;
    !    J, factorization failed on step J, and the solutions could not
    !    be computed.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ) rhs_num

    real ( kind = 8 ) a(n,n+rhs_num)
    real ( kind = 8 ) apivot
    real ( kind = 8 ) factor
    integer ( kind = 4 ) i
    integer ( kind = 4 ) info
    integer ( kind = 4 ) ipivot
    integer ( kind = 4 ) j

    info = 0

    do j = 1, n
       !
       !  Choose a pivot row.
       !
       ipivot = j
       apivot = a(j,j)

       do i = j+1, n
          if ( abs ( apivot ) < abs ( a(i,j) ) ) then
             apivot = a(i,j)
             ipivot = i
          end if
       end do

       if ( apivot == 0.0D+00 ) then
          info = j
          return
       end if
       !
       !  Interchange.
       !
       do i = 1, n + rhs_num
          call r8_swap ( a(ipivot,i), a(j,i) )
       end do
       !
       !  A(J,J) becomes 1.
       !
       a(j,j) = 1.0D+00
       a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
       !
       !  A(I,J) becomes 0.
       !
       do i = 1, n

          if ( i /= j ) then

             factor = a(i,j)
             a(i,j) = 0.0D+00
             a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)

          end if

       end do

    end do

    return
  end subroutine r8mat_solve


  subroutine r8vec_angle_3d ( u, v, angle )

    !*****************************************************************************80
    !
    !! R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) U(3), V(3), the vectors.
    !
    !    Output, real ( kind = 8 ) ANGLE, the angle between the two vectors.
    !
    implicit none

    real ( kind = 8 ) angle
    real ( kind = 8 ) angle_cos
    ! real ( kind = 8 ) r8_acos
    real ( kind = 8 ) u(3)
    real ( kind = 8 ) u_norm
    real ( kind = 8 ) uv_dot
    real ( kind = 8 ) v(3)
    real ( kind = 8 ) v_norm

    uv_dot = dot_product ( u(1:3), v(1:3) )

    u_norm = sqrt ( dot_product ( u(1:3), u(1:3) ) )

    v_norm = sqrt ( dot_product ( v(1:3), v(1:3) ) )

    angle_cos = uv_dot / u_norm / v_norm

    angle = r8_acos ( angle_cos )

    return
  end subroutine r8vec_angle_3d

  subroutine r8vec_cross_3d ( v1, v2, v3 )

    !*****************************************************************************80
    !
    !! R8VEC_CROSS_3D computes the cross product of two vectors in 3D.
    !
    !  Discussion:
    !
    !    The cross product in 3D can be regarded as the determinant of the
    !    symbolic matrix:
    !
    !          |  i  j  k |
    !      det | x1 y1 z1 |
    !          | x2 y2 z2 |
    !
    !      = ( y1 * z2 - z1 * y2 ) * i
    !      + ( z1 * x2 - x1 * z2 ) * j
    !      + ( x1 * y2 - y1 * x2 ) * k
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
    !
    !    Output, real ( kind = 8 ) V3(3), the cross product vector.
    !
    implicit none

    real ( kind = 8 ) v1(3)
    real ( kind = 8 ) v2(3)
    real ( kind = 8 ) v3(3)

    v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
    v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
    v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

    return
  end subroutine r8vec_cross_3d

  function r8vec_length ( dim_num, x )

    !*****************************************************************************80
    !
    !! R8VEC_LENGTH returns the Euclidean length of a vector.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
    !
    !    Input, real ( kind = 8 ) X(DIM_NUM), the vector.
    !
    !    Output, real ( kind = 8 ) R8VEC_LENGTH, the Euclidean length of the vector.
    !
    implicit none

    integer ( kind = 4 ) dim_num

    real ( kind = 8 ) r8vec_length
    real ( kind = 8 ) x(dim_num)

    r8vec_length = sqrt ( sum ( ( x(1:dim_num) )**2 ) )

    return
  end function r8vec_length



  subroutine tetrahedron_centroid ( tetra, centroid )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_CENTROID computes the centroid of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) CENTROID(3), the coordinates of the centroid.
    !
    implicit none

    real ( kind = 8 ) centroid(3)
    integer ( kind = 4 ) i
    real ( kind = 8 ) tetra(3,4)

    do i = 1, 3
       centroid(i) = sum ( tetra(i,1:4) ) / 4.0D+00
    end do

    return
  end subroutine tetrahedron_centroid

  subroutine tetrahedron_circumsphere ( tetra, r, pc )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_CIRCUMSPHERE computes the circumsphere of a tetrahedron.
    !
    !  Discussion:
    !
    !    The circumsphere, or circumscribed sphere, of a tetrahedron is the
    !    sphere that passes through the four vertices.  The circumsphere is
    !    not necessarily the smallest sphere that contains the tetrahedron.
    !
    !    Surprisingly, the diameter of the sphere can be found by solving
    !    a 3 by 3 linear system.  This is because the vectors P2 - P1,
    !    P3 - P1 and P4 - P1 are secants of the sphere, and each forms a
    !    right triangle with the diameter through P1.  Hence, the dot product of
    !    P2 - P1 with that diameter is equal to the square of the length
    !    of P2 - P1, and similarly for P3 - P1 and P4 - P1.  This determines
    !    the diameter vector originating at P1, and hence the radius and
    !    center.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    10 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Adrian Bowyer, John Woodwark,
    !    A Programmer's Geometry,
    !    Butterworths, 1983.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) R, PC(3), the center of the
    !    circumscribed sphere, and its radius.  If the linear system is
    !    singular, then R = -1, PC(1:3) = 0.
    !
    implicit none

    real ( kind = 8 ) a(3,4)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) info
    integer ( kind = 4 ) j
    real ( kind = 8 ) pc(3)
    real ( kind = 8 ) r
    real ( kind = 8 ) tetra(3,4)
    !
    !  Set up the linear system.
    !
    a(1:3,1:3) = transpose ( tetra(1:3,2:4) )

    do j = 1, 3
       a(1:3,j) = a(1:3,j) - tetra(j,1)
    end do

    do i = 1, 3
       a(i,4) = sum ( a(i,1:3)**2 )
    end do
    !
    !  Solve the linear system.
    !
    call r8mat_solve ( 3, 1, a, info )
    !
    !  If the system was singular, return a consolation prize.
    !
    if ( info /= 0 ) then
       r = -1.0D+00
       pc(1:3) = 0.0D+00
       return
    end if
    !
    !  Compute the radius and center.
    !
    r = 0.5D+00 * sqrt ( sum ( a(1:3,4)**2 ) )

    pc(1:3) = tetra(1:3,1) + 0.5D+00 * a(1:3,4)

    return
  end subroutine tetrahedron_circumsphere

  subroutine tetrahedron_dihedral_angles ( tetra, angle )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_DIHEDRAL_ANGLES computes dihedral angles of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) ANGLE(6), the dihedral angles along the
    !    axes AB, AC, AD, BC, BD and CD, respectively.
    !
    implicit none

    real ( kind = 8 ) ab(3)
    real ( kind = 8 ) abc_normal(3)
    real ( kind = 8 ) abd_normal(3)
    real ( kind = 8 ) ac(3)
    real ( kind = 8 ) acd_normal(3)
    real ( kind = 8 ) ad(3)
    real ( kind = 8 ) angle(6)
    real ( kind = 8 ) bc(3)
    real ( kind = 8 ) bcd_normal(3)
    real ( kind = 8 ) bd(3)
    real ( kind = 8 ) cd(3)
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) tetra(3,4)

    call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )

    call r8vec_cross_3d ( ac, ab, abc_normal )
    call r8vec_cross_3d ( ab, ad, abd_normal )
    call r8vec_cross_3d ( ad, ac, acd_normal )
    call r8vec_cross_3d ( bc, bd, bcd_normal )

    call r8vec_angle_3d ( abc_normal, abd_normal, angle(1) )
    call r8vec_angle_3d ( abc_normal, acd_normal, angle(2) )
    call r8vec_angle_3d ( abd_normal, acd_normal, angle(3) )
    call r8vec_angle_3d ( abc_normal, bcd_normal, angle(4) )
    call r8vec_angle_3d ( abd_normal, bcd_normal, angle(5) )
    call r8vec_angle_3d ( acd_normal, bcd_normal, angle(6) )

    angle(1:6) = r8_pi - angle(1:6)

    return
  end subroutine tetrahedron_dihedral_angles

  subroutine tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_EDGES computes the edges of a tetrahedron.
    !
    !  Discussion:
    !
    !    The vertices are A, B, C, D.  The edge from A to B is denoted by AB.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 May 2014
    !
    !  Author:
    !
    !    John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) AB(3), AC(3), AD(3), BC(3), BD(3), CD(3), 
    !    vectors that represent the edges of the tetrahedron.
    !
    implicit none

    real ( kind = 8 ) ab(3)
    real ( kind = 8 ) ac(3)
    real ( kind = 8 ) ad(3)
    real ( kind = 8 ) bc(3)
    real ( kind = 8 ) bd(3)
    real ( kind = 8 ) cd(3)
    real ( kind = 8 ) tetra(3,4)

    ab(1:3) = tetra(1:3,2) - tetra(1:3,1)
    ac(1:3) = tetra(1:3,3) - tetra(1:3,1)
    ad(1:3) = tetra(1:3,4) - tetra(1:3,1)
    bc(1:3) = tetra(1:3,3) - tetra(1:3,2)
    bd(1:3) = tetra(1:3,4) - tetra(1:3,2)
    cd(1:3) = tetra(1:3,4) - tetra(1:3,3)

    return
  end subroutine tetrahedron_edges

  subroutine tetrahedron_edge_length ( tetra, edge_length )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_EDGE_LENGTH returns edge lengths of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) EDGE_LENGTH(6), the length of the edges.
    !
    implicit none

    ! real ( kind = 8 ) r8vec_length
    real ( kind = 8 ) edge_length(6)
    integer ( kind = 4 ) j1
    integer ( kind = 4 ) j2
    integer ( kind = 4 ) k
    real ( kind = 8 ) tetra(3,4)

    k = 0
    do j1 = 1, 3
       do j2 = j1 + 1, 4
          k = k + 1
          edge_length(k) = r8vec_length ( 3, tetra(1:3,j2) - tetra(1:3,j1) )
       end do
    end do

    return
  end subroutine tetrahedron_edge_length

  subroutine tetrahedron_face_angles ( tetra, angles )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_FACE_ANGLES returns the 12 face angles of a tetrahedron.
    !
    !  Discussion:
    !
    !    The tetrahedron has 4 triangular faces.  This routine computes the
    !    3 planar angles associated with each face.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) ANGLES(3,4), the face angles.
    !
    implicit none

    real ( kind = 8 ) angles(3,4)
    real ( kind = 8 ) tri(3,3)
    real ( kind = 8 ) tetra(3,4)
    !
    !  Face 123
    !
    tri(1:3,1:3) = tetra(1:3,1:3)
    call triangle_angles_3d ( tri, angles(1:3,1) )
    !
    !  Face 124
    !
    tri(1:3,1:2) = tetra(1:3,1:2)
    tri(1:3,3) = tetra(1:3,4)
    call triangle_angles_3d ( tri, angles(1:3,2) )
    !
    !  Face 134
    !
    tri(1:3,1) = tetra(1:3,1)
    tri(1:3,2:3) = tetra(1:3,3:4)
    call triangle_angles_3d ( tri, angles(1:3,3) )
    !
    !  Face 234
    !
    tri(1:3,1:3) = tetra(1:3,2:4)
    call triangle_angles_3d ( tri, angles(1:3,4) )

    return
  end subroutine tetrahedron_face_angles

  subroutine tetrahedron_face_areas ( tetra, areas )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_FACE_AREAS returns the 4 face areas of a tetrahedron.
    !
    !  Discussion:
    !
    !    The tetrahedron has 4 triangular faces.  This routine computes the
    !    area of each face.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4) the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) AREAS(4), the face areas.
    !
    implicit none

    real ( kind = 8 ) areas(4)
    real ( kind = 8 ) tri(3,3)
    real ( kind = 8 ) tetra(3,4)
    !
    !  Face 123
    !
    tri(1:3,1:3) = tetra(1:3,1:3)
    call triangle_area_3d ( tri, areas(1) )
    !
    !  Face 124
    !
    tri(1:3,1:2) = tetra(1:3,1:2)
    tri(1:3,3) = tetra(1:3,4)
    call triangle_area_3d ( tri, areas(2) )
    !
    !  Face 134
    !
    tri(1:3,1) = tetra(1:3,1)
    tri(1:3,2:3) = tetra(1:3,3:4)
    call triangle_area_3d ( tri, areas(3) )
    !
    !  Face 234
    !
    tri(1:3,1:3) = tetra(1:3,2:4)
    call triangle_area_3d ( tri, areas(4) )

    return
  end subroutine tetrahedron_face_areas

  subroutine tetrahedron_insphere ( tetra, r, pc )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_INSPHERE finds the insphere of a tetrahedron.
    !
    !  Discussion:
    !
    !    The insphere of a tetrahedron is the inscribed sphere, which touches
    !    each face of the tetrahedron at a single point.
    !
    !    The points of contact are the centroids of the triangular faces
    !    of the tetrahedron.  Therefore, the point of contact for a face
    !    can be computed as the average of the vertices of that face.
    !
    !    The sphere can then be determined as the unique sphere through
    !    the four given centroids.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Philip Schneider, David Eberly,
    !    Geometric Tools for Computer Graphics,
    !    Elsevier, 2002,
    !    ISBN: 1558605940,
    !    LC: T385.G6974.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) R, PC(3), the radius and the center
    !    of the sphere.
    !
    implicit none

    real ( kind = 8 ) b(4,4)
    ! real ( kind = 8 ) r8mat_det_4d
    ! real ( kind = 8 ) r8vec_length
    real ( kind = 8 ) gamma
    real ( kind = 8 ) l123
    real ( kind = 8 ) l124
    real ( kind = 8 ) l134
    real ( kind = 8 ) l234
    real ( kind = 8 ) n123(1:3)
    real ( kind = 8 ) n124(1:3)
    real ( kind = 8 ) n134(1:3)
    real ( kind = 8 ) n234(1:3)
    real ( kind = 8 ) pc(1:3)
    real ( kind = 8 ) r
    real ( kind = 8 ) tetra(1:3,4)
    real ( kind = 8 ) v21(1:3)
    real ( kind = 8 ) v31(1:3)
    real ( kind = 8 ) v41(1:3)
    real ( kind = 8 ) v32(1:3)
    real ( kind = 8 ) v42(1:3)
    real ( kind = 8 ) v43(1:3)

    call tetrahedron_edges ( tetra, v21, v31, v41, v32, v42, v43 )

    call r8vec_cross_3d ( v21, v31, n123 )
    call r8vec_cross_3d ( v41, v21, n124 )
    call r8vec_cross_3d ( v31, v41, n134 )
    call r8vec_cross_3d ( v42, v32, n234 )

    l123 = r8vec_length ( 3, n123 )
    l124 = r8vec_length ( 3, n124 )
    l134 = r8vec_length ( 3, n134 )
    l234 = r8vec_length ( 3, n234 )

    pc(1:3) = ( l234 * tetra(1:3,1)   &
         + l134 * tetra(1:3,2)   &
         + l124 * tetra(1:3,3)   &
         + l123 * tetra(1:3,4) ) &
         / ( l234 + l134 + l124 + l123 )

    b(1:3,1:4) = tetra(1:3,1:4)
    b(4,1:4) = 1.0D+00

    gamma = abs ( r8mat_det_4d ( b ) )

    r = gamma / ( l234 + l134 + l124 + l123 )

    return
  end subroutine tetrahedron_insphere

  subroutine tetrahedron_quality1 ( tetra, quality )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_QUALITY1: "quality" of a tetrahedron.
    !
    !  Discussion:
    !
    !    The quality of a tetrahedron is 3 times the ratio of the radius of
    !    the inscribed sphere divided by that of the circumscribed sphere.
    !
    !    An equilateral tetrahredron achieves the maximum possible quality of 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) QUALITY, the quality of the tetrahedron.
    !
    implicit none

    real ( kind = 8 ) pc(3)
    real ( kind = 8 ) quality
    real ( kind = 8 ) r_in
    real ( kind = 8 ) r_out
    real ( kind = 8 ) tetra(3,4)

    call tetrahedron_circumsphere ( tetra, r_out, pc )

    call tetrahedron_insphere ( tetra, r_in, pc )

    quality = 3.0D+00 * r_in / r_out

    return
  end subroutine tetrahedron_quality1

  subroutine tetrahedron_quality2 ( tetra, quality2 )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_QUALITY2: "quality" of a tetrahedron.
    !
    !  Discussion:
    !
    !    The quality measure #2 of a tetrahedron is:
    !
    !      QUALITY2 = 2 * sqrt ( 6 ) * RIN / LMAX
    !
    !    where
    !
    !      RIN = radius of the inscribed sphere;
    !      LMAX = length of longest side of the tetrahedron.
    !
    !    An equilateral tetrahredron achieves the maximum possible quality of 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    16 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Qiang Du, Desheng Wang,
    !    The Optimal Centroidal Voronoi Tesselations and the Gersho's
    !    Conjecture in the Three-Dimensional Space,
    !    Computers and Mathematics with Applications,
    !    Volume 49, 2005, pages 1355-1373.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the tetrahedron vertices.
    !
    !    Output, real ( kind = 8 ) QUALITY2, the quality of the tetrahedron.
    !
    implicit none

    real ( kind = 8 ) edge_length(6)
    real ( kind = 8 ) l_max
    real ( kind = 8 ) pc(3)
    real ( kind = 8 ) quality2
    real ( kind = 8 ) r_in
    real ( kind = 8 ) tetra(3,4)

    call tetrahedron_edge_length ( tetra, edge_length )

    l_max = maxval ( edge_length(1:6) )

    call tetrahedron_insphere ( tetra, r_in, pc )

    quality2 = 2.0D+00 * sqrt ( 6.0D+00 ) * r_in / l_max

    return
  end subroutine tetrahedron_quality2

  subroutine tetrahedron_quality3 ( tetra, quality3 )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_QUALITY3 computes the mean ratio of a tetrahedron.
    !
    !  Discussion:
    !
    !    This routine computes QUALITY3, the eigenvalue or mean ratio of
    !    a tetrahedron.
    !
    !      QUALITY3 = 12 * ( 3 * volume )^(2/3) / (sum of squares of edge lengths).
    !
    !    This value may be used as a shape quality measure for the tetrahedron.
    !
    !    For an equilateral tetrahedron, the value of this quality measure
    !    will be 1.  For any other tetrahedron, the value will be between
    !    0 and 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 August 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) QUALITY3, the mean ratio of the tetrahedron.
    !
    implicit none

    real ( kind = 8 ) ab(3)
    real ( kind = 8 ) ac(3)
    real ( kind = 8 ) ad(3)
    real ( kind = 8 ) bc(3)
    real ( kind = 8 ) bd(3)
    real ( kind = 8 ) cd(3)
    real ( kind = 8 ) denom
    real ( kind = 8 ) lab
    real ( kind = 8 ) lac
    real ( kind = 8 ) lad
    real ( kind = 8 ) lbc
    real ( kind = 8 ) lbd
    real ( kind = 8 ) lcd
    real ( kind = 8 ) quality3
    real ( kind = 8 ) tetra(3,4)
    real ( kind = 8 ) volume
    !
    !  Compute the vectors representing the sides of the tetrahedron.
    !
    call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )
    !
    !  Compute the squares of the lengths of the sides.
    !
    lab = sum ( ab(1:3)**2 )
    lac = sum ( ac(1:3)**2 )
    lad = sum ( ad(1:3)**2 )
    lbc = sum ( bc(1:3)**2 )
    lbd = sum ( bd(1:3)**2 )
    lcd = sum ( cd(1:3)**2 )
    !
    !  Compute the volume.
    !
    volume = abs ( &
         ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
         + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
         + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

    denom = lab + lac + lad + lbc + lbd + lcd

    if ( denom == 0.0D+00 ) then
       quality3 = 0.0D+00
    else
       quality3 = 12.0D+00 * ( 3.0D+00 * volume )**( 2.0D+00 / 3.0D+00 ) / denom
    end if

    return
  end subroutine tetrahedron_quality3

  subroutine tetrahedron_quality4 ( tetra, quality4 )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_QUALITY4 computes the minimum solid angle of a tetrahedron.
    !
    !  Discussion:
    !
    !    This routine computes a quality measure for a tetrahedron, based
    !    on the sine of half the minimum of the four solid angles.
    !
    !    The quality measure for an equilateral tetrahedron should be 1,
    !    since the solid angles of such a tetrahedron are each equal to pi.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 August 2005
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) QUALITY4, the value of the quality measure.
    !
    implicit none

    real ( kind = 8 ) ab(3)
    real ( kind = 8 ) ac(3)
    real ( kind = 8 ) ad(3)
    real ( kind = 8 ) bc(3)
    real ( kind = 8 ) bd(3)
    real ( kind = 8 ) cd(3)
    real ( kind = 8 ) denom
    real ( kind = 8 ) l1
    real ( kind = 8 ) l2
    real ( kind = 8 ) l3
    real ( kind = 8 ) lab
    real ( kind = 8 ) lac
    real ( kind = 8 ) lad
    real ( kind = 8 ) lbc
    real ( kind = 8 ) lbd
    real ( kind = 8 ) lcd
    real ( kind = 8 ) quality4
    real ( kind = 8 ) tetra(3,4)
    real ( kind = 8 ) volume
    !
    !  Compute the vectors that represent the sides.
    !
    call tetrahedron_edges ( tetra, ab, ac, ad, bc, bd, cd )
    !
    !  Compute the lengths of the sides.
    !
    lab = sqrt ( sum ( ab(1:3)**2 ) )
    lac = sqrt ( sum ( ac(1:3)**2 ) )
    lad = sqrt ( sum ( ad(1:3)**2 ) )
    lbc = sqrt ( sum ( bc(1:3)**2 ) )
    lbd = sqrt ( sum ( bd(1:3)**2 ) )
    lcd = sqrt ( sum ( cd(1:3)**2 ) )
    !
    !  Compute the volume
    !
    volume = abs ( &
         ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
         + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
         + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

    quality4 = 1.0D+00

    l1 = lab + lac
    l2 = lab + lad
    l3 = lac + lad

    denom = ( l1 + lbc ) * ( l1 - lbc ) &
         * ( l2 + lbd ) * ( l2 - lbd ) &
         * ( l3 + lcd ) * ( l3 - lcd )

    if ( denom <= 0.0D+00 ) then
       quality4 = 0.0D+00
    else
       quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
    end if

    l1 = lab + lbc
    l2 = lab + lbd
    l3 = lbc + lbd

    denom = ( l1 + lac ) * ( l1 - lac ) &
         * ( l2 + lad ) * ( l2 - lad ) &
         * ( l3 + lcd ) * ( l3 - lcd )

    if ( denom <= 0.0D+00 ) then
       quality4 = 0.0D+00
    else
       quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
    end if

    l1 = lac + lbc
    l2 = lac + lcd
    l3 = lbc + lcd

    denom = ( l1 + lab ) * ( l1 - lab ) &
         * ( l2 + lad ) * ( l2 - lad ) &
         * ( l3 + lbd ) * ( l3 - lbd )

    if ( denom <= 0.0D+00 ) then
       quality4 = 0.0D+00
    else
       quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
    end if

    l1 = lad + lbd
    l2 = lad + lcd
    l3 = lbd + lcd

    denom = ( l1 + lab ) * ( l1 - lab ) &
         * ( l2 + lac ) * ( l2 - lac ) &
         * ( l3 + lbc ) * ( l3 - lbc )

    if ( denom <= 0.0D+00 ) then
       quality4 = 0.0D+00
    else
       quality4 = min ( quality4, 12.0D+00 * volume / sqrt ( denom ) )
    end if

    quality4 = quality4 * 1.5D+00 * sqrt ( 6.0D+00 )

    return
  end subroutine tetrahedron_quality4

  subroutine tetrahedron_solid_angles ( tetra, angle )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_SOLID_ANGLES computes solid angles of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) ANGLE(4), the solid angles.
    !
    implicit none

    real ( kind = 8 ) angle(4)
    real ( kind = 8 ) dihedral_angles(6)
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) tetra(3,4)

    call tetrahedron_dihedral_angles ( tetra, dihedral_angles )

    angle(1) = dihedral_angles(1) &
         + dihedral_angles(2) &
         + dihedral_angles(3) - r8_pi

    angle(2) = dihedral_angles(1) &
         + dihedral_angles(4) &
         + dihedral_angles(5) - r8_pi

    angle(3) = dihedral_angles(2) &
         + dihedral_angles(4) &
         + dihedral_angles(6) - r8_pi

    angle(4) = dihedral_angles(3) &
         + dihedral_angles(5) &
         + dihedral_angles(6) - r8_pi

    return
  end subroutine tetrahedron_solid_angles

  subroutine tetrahedron_volume ( tetra, volume )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_VOLUME computes the volume of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) VOLUME, the volume of the tetrahedron.
    !
    implicit none

    real ( kind = 8 ) a(4,4)
    ! real ( kind = 8 ) r8mat_det_4d
    real ( kind = 8 ) tetra(3,4)
    real ( kind = 8 ) volume

    a(1:3,1:4) = tetra(1:3,1:4)
    a(4,1:4) = 1.0D+00

    volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

    return
  end subroutine tetrahedron_volume

  subroutine triangle_angles_3d ( t, angle )

    !*****************************************************************************80
    !
    !! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
    !
    !  Discussion:
    !
    !    The law of cosines is used:
    !
    !      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
    !
    !    where GAMMA is the angle opposite side C.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 May 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) ANGLE(3), the angles opposite
    !    sides P1-P2, P2-P3 and P3-P1, in radians.
    !
    implicit none

    real ( kind = 8 ) a
    real ( kind = 8 ) angle(3)
    real ( kind = 8 ) b
    real ( kind = 8 ) c
    ! real ( kind = 8 ) r8_acos
    real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
    real ( kind = 8 ) t(3,3)
    !
    !  Compute the length of each side.
    !
    a = sqrt ( sum ( ( t(1:3,1) - t(1:3,2) )**2 ) )
    b = sqrt ( sum ( ( t(1:3,2) - t(1:3,3) )**2 ) )
    c = sqrt ( sum ( ( t(1:3,3) - t(1:3,1) )**2 ) )
    !
    !  Take care of a ridiculous special case.
    !
    if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
       angle(1:3) = 2.0D+00 * r8_pi / 3.0D+00
       return
    end if

    if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
       angle(1) = r8_pi
    else
       angle(1) = r8_acos ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
    end if

    if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
       angle(2) = r8_pi
    else
       angle(2) = r8_acos ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
    end if

    if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
       angle(3) = r8_pi
    else
       angle(3) = r8_acos ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
    end if

    return
  end subroutine triangle_angles_3d

  subroutine triangle_area_3d ( t, area )

    !*****************************************************************************80
    !
    !! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
    !
    !  Discussion:
    !
    !    This routine uses the fact that the norm of the cross product
    !    of two vectors is the area of the parallelogram they form.
    !
    !    Therefore, the area of the triangle is half of that value.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Adrian Bowyer, John Woodwark,
    !    A Programmer's Geometry,
    !    Butterworths, 1983.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) AREA, the area of the triangle.
    !
    implicit none

    real ( kind = 8 ) area
    real ( kind = 8 ) cross(3)
    real ( kind = 8 ) t(3,3)
    !
    !  Compute the cross product vector.
    !
    cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
         - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

    cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
         - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

    cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
         - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

    area = 0.5D+00 * sqrt ( sum ( cross(1:3)**2 ) )

    return
  end subroutine triangle_area_3d

end module tet_props


! program main
!   use tet_props
!   implicit none

!   print *, 'Hello!'

!   ! done here
! end program main

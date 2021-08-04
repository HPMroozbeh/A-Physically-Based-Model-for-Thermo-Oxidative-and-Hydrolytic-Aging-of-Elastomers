function [ D , W ] = integration_points_data(Ndir)
if     Ndir == 8
    a1 = 1/sqrt(3) ;
    a2 = 1/sqrt(3) ;
    a3 = 1/sqrt(3) ;
    w1 = 1/8 ;
    w2 = 1/8 ;
    w3 = 1/8 ;
    point_information = ...
        [ a1   a2   a3  w1
        a1  -a2   a3  w2
        -a1   a2   a3  w3
        -a1  -a2   a3  w3];
    
    D = [point_information(1:4,1:3) ; -point_information(1:4,1:3)];
    W = [point_information(1:4,4)   ; point_information(1:4,4)];
elseif Ndir == 42
    a1 = 0.707106781187 ;
    a2 = 0.387907304067 ;
    a3 = 0.836095596749 ;
    w1 = 0.0265214244093 ;
    w2 = 0.0199301476312 ;
    w3 = 0.0250712367487 ;
    
    point_information = ...
        [ 1  0  0  w1
        0  1  0  w1
        0  0  1  w1
        a1  a1  0  w2
        a1 -a1  0  w2
        a1  0  a1  w2
        a1  0 -a1  w2
        0  a1  a1  w2
        0  a1 -a1  w2
        a2  a2  a3  w3
        a2  a2 -a3  w3
        a2 -a2  a3  w3
        a2 -a2 -a3  w3
        a2  a3  a2  w3
        a2  a3 -a2  w3
        a2 -a3  a2  w3
        a2 -a3 -a2  w3
        a3  a2  a2  w3
        a3  a2 -a2  w3
        a3 -a2  a2  w3
        a3 -a2 -a2  w3] ;
    
    D = [point_information(1:21,1:3) ; -point_information(1:21,1:3)];
    W = [point_information(1:21,4)   ; point_information(1:21,4)];
end
end
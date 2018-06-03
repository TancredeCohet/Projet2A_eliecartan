%> @file  build_Q2Pressure.m
%> @filebrief Extrapolates pressure on Q2 mesh
%======================================================================
%> @brief Extrapolates the pressure on the Q2 mesh so that it can be easily loaded and used into Paraview
%> @param t : the connectivity matrix of the domain's mesh (Q1)
%> @param P : the pressure described on the Q1 mesh
%> @param nb_dof : the number of Q2 points
%> @retval P2 : the pressure described on the Q2 mesh
function [P2] = build_Q2Pressure(t,P,nb_dof)
%***********************************************************************
% Copyright ©2013 Inria - Université de Lorraine
% This software is under GNU General Public Licence
% (http://www.gnu.org/licenses/gpl-3.0.txt)
%***********************************************************************
% Construit la pression sur le maillage Q2 pour la visualisation VTK
% t : la matrice de connectivit?? de la grille Q1
% P : la pression d??crite sur le maillage Q1
% variable de retour P2 : la pression d??crite sur le maillage Q2
P2 = zeros(nb_dof,1);

P2(t(:,1:8))=P(t(:,1:8));
P2(t(:,9))=(P(t(:,1))+P(t(:,2)))/2;
P2(t(:,10))=(P(t(:,2))+P(t(:,3)))/2;
P2(t(:,11))=(P(t(:,3))+P(t(:,4)))/2;
P2(t(:,12))=(P(t(:,1))+P(t(:,4)))/2;

P2(t(:,13))=(P(t(:,1))+P(t(:,5)))/2;
P2(t(:,14))=(P(t(:,2))+P(t(:,6)))/2;
P2(t(:,15))=(P(t(:,3))+P(t(:,7)))/2;
P2(t(:,16))=(P(t(:,4))+P(t(:,8)))/2;

P2(t(:,17))=(P(t(:,5))+P(t(:,6)))/2;
P2(t(:,18))=(P(t(:,6))+P(t(:,7)))/2;
P2(t(:,19))=(P(t(:,7))+P(t(:,8)))/2;
P2(t(:,20))=(P(t(:,5))+P(t(:,8)))/2;

P2(t(:,21))=(P(t(:,1))+P(t(:,2))+P(t(:,3))+P(t(:,4)))/4;
P2(t(:,22))=(P(t(:,1))+P(t(:,2))+P(t(:,5))+P(t(:,6)))/4;
P2(t(:,23))=(P(t(:,2))+P(t(:,3))+P(t(:,6))+P(t(:,7)))/4;
P2(t(:,24))=(P(t(:,3))+P(t(:,4))+P(t(:,7))+P(t(:,8)))/4;
P2(t(:,25))=(P(t(:,1))+P(t(:,4))+P(t(:,5))+P(t(:,8)))/4;
P2(t(:,26))=(P(t(:,5))+P(t(:,6))+P(t(:,7))+P(t(:,8)))/4;

P2(t(:,27))=(P(t(:,1))+P(t(:,2))+P(t(:,3))+P(t(:,4))+P(t(:,5))+P(t(:,6))+P(t(:,7))+P(t(:,8)))/8;


end


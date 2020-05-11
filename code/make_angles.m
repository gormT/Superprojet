function angles = make_angles(angles_triangles,angle_steps,Nbtri,N)
% Donne tous les angles utilisés lors de l'intégration numérique.

angles = zeros(Nbtri,(N+1)*3);

angles(:,[1 N+1 N+2 2*N+2 2*N+3 3*N+3]) = angles_triangles;
for t=1:N-1
    angles(:,[t+1 t+N+2 t+2*N+3]) = angles(:,[t t+N+1 t+2*N+2]).*angle_steps;
end

end
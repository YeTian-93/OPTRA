function newLap = AccGossip(Lap, K, c)
global Num_Nodes
a = 1;   new_a = c;
Z = eye(Num_Nodes);   new_Z = c*(eye(Num_Nodes) - Lap);
for k = 1:K-1
    temp_a = 2*c*new_a - a;
    temp_Z = 2*c*(eye(Num_Nodes) - Lap)*new_Z - Z;
    a = new_a; Z = new_Z;
    new_a = temp_a; new_Z = temp_Z;
end
newLap = eye(Num_Nodes) - new_Z/new_a;
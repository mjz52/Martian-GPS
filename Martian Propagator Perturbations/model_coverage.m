% Michael Zakoworotny
% Testing circle of coverage of a satellite on surface

function model_coverage(r, geod, alpha)

x = r(:,1); y = r(:,2); z = r(:,3);
lon = geod(:,1); lat = geod(:,2); h = (:,3);
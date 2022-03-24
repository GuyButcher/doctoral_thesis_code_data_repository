function extension = force_domain_unfolding(forces,thresholds)
    for i = 1:size(forces,2)
        extension(i) = subDomainSpring(forces(i),thresholds(1));
        for j = 2:size(thresholds,2)
            extension(i) = extension(i) + subDomainSpring(forces(i),thresholds(j));
        end
    end
end


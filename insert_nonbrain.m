function new_cifti = insert_nonbrain(input, hem, template)
count = 1;
new_cifti = [];
if strcmp(hem, 'both')
    for i = 1:length(template.brainstructure)
        if template.brainstructure(i)>0
            new_cifti(i,:) = input(count,:);
            count = count + 1;
        else
            new_cifti(i,:) = 0;
        end
    end
else
    left_hem = template.brainstructure(1:32492);
    right_hem = template.brainstructure(32493:end);
    if strcmp(hem, 'left')
        for l = 1:length(left_hem)
            if left_hem(l)>0
                new_cifti(l,:) = input(count, :);
                count = count + 1;
            else
                new_cifti(l,:) = 0;
            end
        end
    elseif strcmp(hem, 'right')
        for r = 1:length(right_hem)
            if right_hem(r)>0
                new_cifti(r,:) = input(count,:);
                count = count + 1;
            else
                new_cifti(r,:) = 0;
            end
        end
    end
end
end
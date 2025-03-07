
figure(1000)
clf
tiledlayout(3,3)

for orientation = 1:3
    for comp = 1:3
        nexttile
        hold on
        for ii = 1:length(models)
            plot(models(ii).coil_detectors(orientation).coil_system.coils{1, 1}.B_values(:,comp),'DisplayName',sprintf('m %d, or %d, c %d', ii, orientation, comp))
        end
        legend
    end
end

figure(2000)
clf
tiledlayout(3,1)
for orientation = 1:3

        nexttile
        hold on
        for ii = 1:length(models)
            plot(models(ii).coil_detectors(orientation).coil_system.coils{1, 1}.B_values(:,comp),'DisplayName',sprintf('m %d, or %d, c %d', ii, orientation, comp))
        end
        legend

end


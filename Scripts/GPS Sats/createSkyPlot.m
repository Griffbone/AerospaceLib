function createSkyPlot(azimuths, elevations, labels)
    % Convert azimuths and elevations to polar coordinates
    thetas = deg2rad(90 - azimuths);
    rs = 90 - elevations;

    % Create polar plot
    polarplot(thetas, rs, 'o'); % Plot points

    % Add labels at each point
    for i = 1:numel(azimuths)
        text(thetas(i), rs(i), labels{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end

    % Add cardinal direction labels
    cardinal_labels = {'N', 'E', 'S', 'W'};
    for i = 1:4
        theta_cardinal = deg2rad((i - 1) * 90); % Angle for cardinal direction
        polarplot([theta_cardinal, theta_cardinal], [0, 90], '-k'); % Plot line
        text(theta_cardinal, 95, cardinal_labels{i}, 'HorizontalAlignment', 'center');
    end

    % Set polar plot properties
    set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');
    title('Skyplot');
end

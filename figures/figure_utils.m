% FIGURE_UTILS  The colors, type and line defaults, panel letters and export that the figure
% scripts share.
%
% Implements : the style and export of Figs. 1 to 6, S1 and S2
% Inputs     : none (a class of constant colors and static methods)
% Outputs    : colors (Okabe-Ito and the project's), the type constants FONT, SIZE_SMALL,
%              SIZE_LABEL and SIZE_LETTER, setup_defaults (Arial 9 pt, 1.5 pt lines, ticks out, no
%              box, for every new figure), create_figure (a figure of a given size in inches with
%              the same defaults), add_panel_label, style_axis, export (TIFF, PDF and PNG, after
%              validate_figure), architecture_palette, style_panel (0.6 pt axes, 7 and 8 pt type),
%              head_cm (a panel letter and title), pass_fail, check_pdf_width, export_tiff (LZW)
% Seed       : none (no random draws)
%
% Usage: figure_utils.setup_defaults(); fig = figure_utils.create_figure(7.5, 5.0);
% figure_utils.add_panel_label(gca, 'A'); figure_utils.export(fig, 'Fig1').

classdef figure_utils

    properties (Constant)
        % Okabe-Ito colorblind-safe palette
        BLUE        = [0.000, 0.447, 0.698];
        ORANGE      = [0.902, 0.624, 0.000];
        GREEN       = [0.000, 0.620, 0.451];
        VERMILLION  = [0.835, 0.369, 0.000];
        SKY_BLUE    = [0.337, 0.706, 0.914];
        PURPLE      = [0.800, 0.475, 0.655];
        BLACK       = [0.000, 0.000, 0.000];
        GRAY        = [0.500, 0.500, 0.500];
        INK         = [0.150, 0.150, 0.150];   % observed and control arms; axes and text
        NOTE        = [0.350, 0.350, 0.350];   % notes and secondary marks
        DIAGONAL_GRAY = [0.502, 0.502, 0.502]; % Fig. 6A's diagonal
        % The game-state fills: Fig. 1C's "dilemma holds" / "dilemma escaped" cells and Fig. 2A's selection strip
        DILEMMA_HOLDS_TINT   = [0.878, 0.749, 0.706];   % #E0BFB4
        DILEMMA_ESCAPED_TINT = [0.757, 0.808, 0.914];   % #C1CEE9

        % Semantic color assignments for this project (4-gene model)
        COLOR_COOP     = [0.000, 0.447, 0.698];  % Blue
        COLOR_DEFECT   = [0.835, 0.369, 0.000];  % Vermillion

        % Full palette as matrix (for color order)
        PALETTE = [
            0.000, 0.447, 0.698;   % Blue
            0.902, 0.624, 0.000;   % Orange
            0.000, 0.620, 0.451;   % Bluish Green
            0.835, 0.369, 0.000;   % Vermillion
            0.337, 0.706, 0.914;   % Sky Blue
            0.800, 0.475, 0.655;   % Reddish Purple
            0.941, 0.894, 0.259;   % Yellow
        ];

        % Type: one font; small text (ticks, legends, direct labels), axis and panel titles,
        % panel letters (MATLAB FontSize)
        FONT         = 'Arial';
        SIZE_SMALL   = 7;
        SIZE_LABEL   = 8;
        SIZE_LETTER  = 10;
    end

    methods (Static)

        function setup_defaults()
            %% Set global defaults for all subsequent figures
            set(groot, 'DefaultAxesFontName', figure_utils.FONT);
            set(groot, 'DefaultAxesFontSize', 9);
            set(groot, 'DefaultTextFontName', figure_utils.FONT);
            set(groot, 'DefaultTextFontSize', 9);
            set(groot, 'DefaultLineLineWidth', 1.5);
            set(groot, 'DefaultLineMarkerSize', 6);
            set(groot, 'DefaultAxesLineWidth', 0.75);
            set(groot, 'DefaultAxesBox', 'off');
            set(groot, 'DefaultAxesTickDir', 'out');
            set(groot, 'DefaultAxesTickLength', [0.015, 0.015]);
            set(groot, 'DefaultAxesXColor', figure_utils.INK);
            set(groot, 'DefaultAxesYColor', figure_utils.INK);
            set(groot, 'DefaultAxesColorOrder', figure_utils.PALETTE);
            set(groot, 'DefaultFigureColor', 'w');
            set(groot, 'DefaultAxesColor', 'w');
        end

        function fig = create_figure(width_inches, height_inches)
            %% Create a publication-ready figure with correct dimensions
            if nargin < 1, width_inches = 7.5; end
            if nargin < 2, height_inches = 5.0; end

            fig = figure('Color', 'w');
            fig.Units = 'inches';
            fig.Position = [1, 1, width_inches, height_inches];

            % Set figure-level defaults
            set(fig, 'DefaultAxesFontName', figure_utils.FONT);
            set(fig, 'DefaultAxesFontSize', 9);
            set(fig, 'DefaultTextFontName', figure_utils.FONT);
            set(fig, 'DefaultTextFontSize', 9);
            set(fig, 'DefaultLineLineWidth', 1.5);
            set(fig, 'DefaultLineMarkerSize', 6);
            set(fig, 'DefaultAxesLineWidth', 0.75);
            set(fig, 'DefaultAxesBox', 'off');
            set(fig, 'DefaultAxesTickDir', 'out');
            set(fig, 'DefaultAxesTickLength', [0.015, 0.015]);
            set(fig, 'DefaultAxesXColor', figure_utils.INK);
            set(fig, 'DefaultAxesYColor', figure_utils.INK);
            set(fig, 'DefaultAxesColorOrder', figure_utils.PALETTE);
        end

        function add_panel_label(ax, label, pos)
            %% Add bold panel label (A, B, C...) in upper-left corner.
            %  Default position (-0.10, 1.03) sits just outside the axes;
            %  pass pos=[x y] (normalized) to override (e.g. inside axes).
            if nargin < 3 || isempty(pos), pos = [-0.10, 1.03]; end
            text(ax, pos(1), pos(2), label, ...
                'Units', 'normalized', ...
                'FontSize', 12, ...
                'FontWeight', 'bold', ...
                'FontName', figure_utils.FONT, ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom');
        end

        function style_axis(ax)
            %% Apply consistent styling to an axis
            set(ax, 'FontName', figure_utils.FONT, 'FontSize', 9);
            set(ax, 'TickDir', 'out', 'TickLength', [0.015, 0.015]);
            set(ax, 'Box', 'off');
            set(ax, 'LineWidth', 0.75);
            set(ax, 'XColor', figure_utils.INK);
            set(ax, 'YColor', figure_utils.INK);
        end

        function export(fig, filename, dpi, spec, padding)
            %% Export a figure in three formats:
            %    <filename>.tif   LZW-compressed TIFF at dpi (default 300)
            %    <filename>.pdf   vector PDF
            %    <filename>.png   raster at 150 dpi, a preview
            %  It first runs validate_figure(fig, spec) and warns on any issue, then
            %  fixes the axis limits and ticks and sets the paper size to the figure.
            %  padding (default 'tight') is exportgraphics's Padding: 'tight' crops to
            %  the drawn content, a number (points) adds that much space on every side.
            %  After export it warns when the TIFF is over 10 MB, narrower than 789 or
            %  wider than 2250 pixels, taller than 2625 pixels, not RGB truecolor, not
            %  24-bit, or has an alpha channel.
            %  Known limitation, not handled here: MATLAB sets the vector PDF's text at
            %  whole-pixel sizes and spacing, so the TIFF's text can run slightly
            %  larger than the PDF's; check the TIFF.
            if nargin < 3, dpi = 300; end
            if nargin < 4, spec = struct(); end
            if nargin < 5 || isempty(padding), padding = 'tight'; end

            % Auto-populate spec.name from filename
            [~, basename] = fileparts(filename);
            if ~isfield(spec, 'name'), spec.name = basename; end

            % Run validation before export
            [valOk, valIssues] = validate_figure(fig, spec);
            if ~valOk
                warning('figure_utils:issues', ...
                    '%s has %d validation issue(s) -- exporting anyway.', ...
                    basename, numel(valIssues));
            end

            % Freeze axis limits
            all_axes = findall(fig, 'Type', 'axes');
            for ax = all_axes'
                set(ax, 'XLimMode', 'manual', 'YLimMode', 'manual');
                set(ax, 'XTickMode', 'manual', 'YTickMode', 'manual');
            end

            % Set paper size to match figure size
            fig.PaperPositionMode = 'auto';
            fig_pos = fig.Position;
            fig.PaperUnits = 'inches';
            fig.PaperSize = [fig_pos(3), fig_pos(4)];
            fig.PaperPosition = [0, 0, fig_pos(3), fig_pos(4)];

            % Export TIFF with LZW compression: exportgraphics writes no compressed
            % TIFF, so the figure is rendered to an intermediate PNG and re-encoded
            % with imwrite and LZW
            tiff_file = [filename, '.tif'];
            tmp_png = [tempname, '.png'];
            try
                exportgraphics(fig, tmp_png, 'Resolution', dpi, 'ContentType', 'image', 'Padding', padding);
                img = imread(tmp_png);
                imwrite(img, tiff_file, 'tif', ...
                    'Compression', 'lzw', 'Resolution', dpi);
                delete(tmp_png);
            catch e
                if exist(tmp_png, 'file'), delete(tmp_png); end
                warning('TIFF export failed: %s. Falling back to exportgraphics.', e.message);
                exportgraphics(fig, tiff_file, 'Resolution', dpi, 'ContentType', 'image', 'Padding', padding);
            end

            % Export PDF (for review/proofing)
            pdf_file = [filename, '.pdf'];
            w = warning('off', 'MATLAB:print:ContentTypeImageSuggested');
            exportgraphics(fig, pdf_file, 'ContentType', 'vector', 'Padding', padding);
            warning(w);

            % PNG preview at 150 dpi.
            png_file = [filename, '.png'];
            exportgraphics(fig, png_file, 'Resolution', 150, 'Padding', padding);

            % Verify file size and TIFF dimensions
            info = dir(tiff_file);
            size_mb = info.bytes / 1e6;
            if size_mb > 10
                warning('TIFF %.1f MB exceeds the 10 MB limit', size_mb);
            end

            % Post-export TIFF verification (dimensions, color, metadata)
            try
                tiffInfo = imfinfo(tiff_file);
                pixW = tiffInfo.Width;
                pixH = tiffInfo.Height;
                if pixW < 789
                    warning('TIFF width %d px < min 789 px', pixW);
                end
                if pixW > 2250
                    warning('TIFF width %d px > max 2250 px', pixW);
                end
                if pixH > 2625
                    warning('TIFF height %d px > max 2625 px', pixH);
                end
                % RGB 8-bit, no alpha channel
                if ~strcmp(tiffInfo.ColorType, 'truecolor')
                    warning('TIFF ColorType is "%s"; RGB truecolor required', tiffInfo.ColorType);
                end
                if tiffInfo.BitDepth ~= 24
                    warning('TIFF BitDepth is %d; 24-bit (8-bit RGB) required', tiffInfo.BitDepth);
                end
                if isfield(tiffInfo, 'HasAlpha') && tiffInfo.HasAlpha
                    warning('TIFF has an alpha channel; none allowed');
                end
            catch
                % imfinfo may fail on some TIFF variants; skip silently
            end

            fprintf('Exported: %s (%.2f MB) + %s + %s\n', tiff_file, size_mb, pdf_file, png_file);
        end

        function p = architecture_palette()
            %% The signaling-architecture palette, in the order graded separate,
            %  graded contrast, saturating separate, saturating contrast: purple,
            %  orange, green, ink. Blue and vermillion are kept for cooperation and
            %  defection, so no architecture wears them. Returns .colors (4x3),
            %  .linestyles (dashed separate, solid contrast), .markers (square
            %  graded, circle saturating), .labels (unhyphenated, as the paper names
            %  them), .short_labels, one struct per architecture
            %  (.graded_separate, ...) with .color, .linestyle and .marker, and
            %  .spatial (vermillion), which no figure draws.
            p.colors = [
                figure_utils.PURPLE;        % graded-separate
                figure_utils.ORANGE;        % graded-contrast
                figure_utils.GREEN;         % saturating-separate
                figure_utils.INK;           % saturating-contrast
            ];
            p.linestyles = {'--', '-', '--', '-'};
            p.markers    = {'s', 's', 'o', 'o'};   % graded=square, saturating=circle
            p.labels     = {'graded separate', 'graded contrast', ...
                            'saturating separate', 'saturating contrast'};
            p.short_labels = {'grad-sep', 'grad-con', 'sat-sep', 'sat-con'};

            p.graded_separate     = struct('color', p.colors(1,:), ...
                'linestyle', p.linestyles{1}, 'marker', p.markers{1});
            p.graded_contrast     = struct('color', p.colors(2,:), ...
                'linestyle', p.linestyles{2}, 'marker', p.markers{2});
            p.saturating_separate = struct('color', p.colors(3,:), ...
                'linestyle', p.linestyles{3}, 'marker', p.markers{3});
            p.saturating_contrast = struct('color', p.colors(4,:), ...
                'linestyle', p.linestyles{4}, 'marker', p.markers{4});

            % No figure draws it.
            p.spatial = struct('color', figure_utils.VERMILLION, ...
                'linestyle', '-', 'marker', 's');
        end

        function style_panel(ax)
            %% 0.6 pt axes in ink, ticks out, no box, no grid; SIZE_SMALL ticks and SIZE_LABEL
            %  axis titles (the label multiplier pinned to 1). Figs. 6 and S2.
            set(ax, 'FontName', figure_utils.FONT, 'FontSize', figure_utils.SIZE_SMALL, 'LineWidth', 0.6, ...
                'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off', 'XColor', figure_utils.INK, ...
                'YColor', figure_utils.INK, 'XTickLabelRotation', 0, 'YTickLabelRotation', 0, ...
                'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
            grid(ax, 'off');
            set([ax.XLabel ax.YLabel], 'FontSize', figure_utils.SIZE_LABEL, 'Color', figure_utils.INK, ...
                'FontName', figure_utils.FONT);
        end

        function head_cm(ax, letter, ttl, dxCm)
            %% Panel letter (SIZE_LETTER, bold) dxCm left of the plot box and the panel title
            %  (SIZE_LABEL) on one baseline above it. Used by Fig. 6.
            u = ax.Units; ax.Units = 'centimeters'; wcm = ax.Position(3); ax.Units = u;
            text(ax, -dxCm / wcm, 1.04, letter, 'Units', 'normalized', 'FontSize', figure_utils.SIZE_LETTER, ...
                'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            text(ax, 0, 1.04, ttl, 'Units', 'normalized', 'FontSize', figure_utils.SIZE_LABEL, ...
                'Color', figure_utils.INK, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        end

        function s = pass_fail(t)
            %% 'PASS' or 'FAIL', for the printed checks.
            if t, s = 'PASS'; else, s = 'FAIL'; end
        end

        function check_pdf_width(pdfFile, widthCm)
            %% Stops unless the exported PDF's MediaBox width equals widthCm within 2 pt; prints it.
            txt = fileread(pdfFile);
            tok = regexp(txt, '/MediaBox\s*\[\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\]', 'tokens', 'once');
            if isempty(tok), error('figure_utils:pdf_width', 'no MediaBox in %s', pdfFile); end
            v = str2double(tok); wPt = v(3) - v(1); target = widthCm / 2.54 * 72;
            fprintf('  CONTROL PDF MediaBox width %.1f pt (target %.1f, tolerance 2): %s\n', wPt, target, ...
                figure_utils.pass_fail(abs(wPt - target) <= 2));
            if abs(wPt - target) > 2
                error('figure_utils:pdf_width', 'PDF width %.1f pt differs from %.1f pt', wPt, target);
            end
        end

        function export_tiff(fig, tifPath, dpi)
            %% RGB TIFF with LZW compression: exportgraphics writes no compressed TIFF, so the
            %  figure is rendered to a temporary PNG at dpi (default 600; deleted afterwards, also
            %  on error) and re-encoded by imwrite. Figs. 1 and 3.
            if nargin < 3 || isempty(dpi), dpi = 600; end
            tmp_png = [tempname(), '.png'];
            try
                exportgraphics(fig, tmp_png, 'Resolution', dpi, 'ContentType', 'image');
                img = imread(tmp_png);
                imwrite(img, tifPath, 'tif', 'Compression', 'lzw', 'Resolution', dpi);
                delete(tmp_png);
            catch e
                if exist(tmp_png, 'file'), delete(tmp_png); end
                rethrow(e);
            end
        end

    end
end

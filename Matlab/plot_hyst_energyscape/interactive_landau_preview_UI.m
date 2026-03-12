function interactive_landau_preview_UI()
% INTERACTIVE_LANDAU_PREVIEW_UI
% Slider-based UI to tune (a1,b1,a2,b2,c2) and preview:
%   Fpsi(psi) = 0.5*a1*psi^2 + (1/4!)*b1*psi^4
%   FX(X)     = 0.5*a2*X^2  + (1/4!)*b2*X^4 + (1/6!)*c2*X^6

    % ===================== defaults =====================
    p = struct();
    p.a1 = 0.2;   p.b1 = 1.0;
    p.a2 = 0.5;   p.b2 = -3.0;  p.c2 = 1.0;

    view = struct();
    view.psi_max  = 20.0;
    view.X_max    = 10.0;
    view.Ngrid    = 401;
    view.shiftMin = true;

    % slider ranges (you can edit as you like)
    r = struct();
    r.a1 = [-50, 50];
    r.b1 = [  0, 20];     % usually positive for stability
    r.a2 = [-50, 50];
    r.b2 = [-50, 50];     % can be negative to create double well
    r.c2 = [  0, 20];     % usually positive

    % ===================== UI layout =====================
    fig = uifigure('Name','Landau Preview (a1,b1,a2,b2,c2)','Position',[100 80 1100 620]);

    gl = uigridlayout(fig,[1 2]);
    gl.ColumnWidth = {340, '1x'};
    gl.RowHeight   = {'1x'};

    % ---- left panel: controls ----
    pnl = uipanel(gl,'Title','Parameters','FontWeight','bold');
    gL  = uigridlayout(pnl,[14 3]);
    gL.ColumnWidth = {70,'1x',70};
    gL.RowHeight   = repmat({28},1,14);
    gL.Padding = [10 10 10 10];
    gL.RowSpacing = 6;
    gL.ColumnSpacing = 8;

    % helper to add a slider row
    [sldA1, edtA1] = add_slider_row_(gL, 1, 'a1', r.a1, p.a1);
    [sldB1, edtB1] = add_slider_row_(gL, 2, 'b1', r.b1, p.b1);
    [sldA2, edtA2] = add_slider_row_(gL, 3, 'a2', r.a2, p.a2);
    [sldB2, edtB2] = add_slider_row_(gL, 4, 'b2', r.b2, p.b2);
    [sldC2, edtC2] = add_slider_row_(gL, 5, 'c2', r.c2, p.c2);

    % view settings
    uilabel(gL,'Text','psiMax','HorizontalAlignment','right'); 
    edtPsiMax = uieditfield(gL,'numeric','Value',view.psi_max,'Limits',[0.1, 1e6]);
    uilabel(gL,'Text',''); %#ok<*NASGU>

    uilabel(gL,'Text','XMax','HorizontalAlignment','right'); 
    edtXMax = uieditfield(gL,'numeric','Value',view.X_max,'Limits',[0.1, 1e6]);
    uilabel(gL,'Text','');

    uilabel(gL,'Text','Ngrid','HorizontalAlignment','right'); 
    edtNgrid = uieditfield(gL,'numeric','Value',view.Ngrid,'Limits',[101, 20001],'RoundFractionalValues','on');
    uilabel(gL,'Text','');

    cbShift = uicheckbox(gL,'Text','shift by min (F-min)','Value',view.shiftMin);
    cbShift.Layout.Row = 9; cbShift.Layout.Column = [1 3];

    btnReset = uibutton(gL,'Text','Reset defaults');
    btnReset.Layout.Row = 10; btnReset.Layout.Column = [1 3];

    txtInfo = uitextarea(gL,'Editable','off','Value',{'Ready.'});
    txtInfo.Layout.Row = [11 14];
    txtInfo.Layout.Column = [1 3];

    % ---- right panel: plots ----
    pnlP = uipanel(gl,'Title','Preview','FontWeight','bold');
    gR = uigridlayout(pnlP,[1 1]);
    gR.Padding = [10 10 10 10];

    % use tiledlayout inside a uipanel with standard axes
    axContainer = uipanel(gR,'BorderType','none');
    axContainer.Layout.Row = 1; axContainer.Layout.Column = 1;

    % create normal figure-like axes embedded in uifigure
    % (uiaxes supports basic plot well; latex works too)
    t = tiledlayout(axContainer,1,2,'TileSpacing','compact','Padding','compact');

    ax1 = nexttile(t,1);
    ax2 = nexttile(t,2);
    ax1.TickLabelInterpreter = 'latex';
    ax2.TickLabelInterpreter = 'latex';

    % initialize line objects
    [psi, Fpsi_v, X, FX_v] = eval_curves_(p, view);
    h1 = plot(ax1, psi, Fpsi_v, 'LineWidth', 2); grid(ax1,'on'); box(ax1,'on');
    h2 = plot(ax2, X,   FX_v,   'LineWidth', 2); grid(ax2,'on'); box(ax2,'on');

    xlabel(ax1,'$\psi$','Interpreter','latex');
    ylabel(ax1, view_ylabel_(view.shiftMin, "psi"), 'Interpreter','latex');
    title(ax1, title_psi_(p), 'Interpreter','latex','FontWeight','normal');

    xlabel(ax2,'$X$','Interpreter','latex');
    ylabel(ax2, view_ylabel_(view.shiftMin, "X"), 'Interpreter','latex');
    title(ax2, title_X_(p), 'Interpreter','latex','FontWeight','normal');

    % ===================== callbacks wiring =====================
    % sliders: update live while dragging + after release
    hook_slider_(sldA1, edtA1, @(v)set_param_("a1",v));
    hook_slider_(sldB1, edtB1, @(v)set_param_("b1",v));
    hook_slider_(sldA2, edtA2, @(v)set_param_("a2",v));
    hook_slider_(sldB2, edtB2, @(v)set_param_("b2",v));
    hook_slider_(sldC2, edtC2, @(v)set_param_("c2",v));

    % edit fields manual entry: also update slider + plot
    hook_edit_(edtA1, sldA1, @(v)set_param_("a1",v));
    hook_edit_(edtB1, sldB1, @(v)set_param_("b1",v));
    hook_edit_(edtA2, sldA2, @(v)set_param_("a2",v));
    hook_edit_(edtB2, sldB2, @(v)set_param_("b2",v));
    hook_edit_(edtC2, sldC2, @(v)set_param_("c2",v));

    % view settings changes
    edtPsiMax.ValueChangedFcn = @(~,~)set_view_();
    edtXMax.ValueChangedFcn   = @(~,~)set_view_();
    edtNgrid.ValueChangedFcn  = @(~,~)set_view_();
    cbShift.ValueChangedFcn   = @(~,~)set_view_();

    btnReset.ButtonPushedFcn  = @(~,~)do_reset_();

    % ===================== nested helper funcs =====================
    function set_param_(name, v)
        p.(name) = v;
        refresh_plot_(false);
    end

    function set_view_()
        view.psi_max  = edtPsiMax.Value;
        view.X_max    = edtXMax.Value;
        view.Ngrid    = round(edtNgrid.Value);
        view.shiftMin = logical(cbShift.Value);
        refresh_plot_(true);
    end

    function do_reset_()
        p.a1 = 0.2; p.b1 = 1.0;
        p.a2 = 0.5; p.b2 = -3.0; p.c2 = 1.0;

        view.psi_max  = 10.0;
        view.X_max    = 10.0;
        view.Ngrid    = 2001;
        view.shiftMin = true;

        % push values back to UI
        set_slider_edit_(sldA1, edtA1, p.a1);
        set_slider_edit_(sldB1, edtB1, p.b1);
        set_slider_edit_(sldA2, edtA2, p.a2);
        set_slider_edit_(sldB2, edtB2, p.b2);
        set_slider_edit_(sldC2, edtC2, p.c2);

        edtPsiMax.Value = view.psi_max;
        edtXMax.Value   = view.X_max;
        edtNgrid.Value  = view.Ngrid;
        cbShift.Value   = view.shiftMin;

        refresh_plot_(true);
    end

    function refresh_plot_(update_axes_labels)
        try
            [psi_, Fpsi_, X_, FX_] = eval_curves_(p, view);

            % update line data (fast, smooth)
            h1.XData = psi_;
            h1.YData = Fpsi_;
            h2.XData = X_;
            h2.YData = FX_;

            if update_axes_labels
                ylabel(ax1, view_ylabel_(view.shiftMin,"psi"), 'Interpreter','latex');
                ylabel(ax2, view_ylabel_(view.shiftMin,"X"),   'Interpreter','latex');
            end

            title(ax1, title_psi_(p), 'Interpreter','latex','FontWeight','normal');
            title(ax2, title_X_(p),   'Interpreter','latex','FontWeight','normal');

            txtInfo.Value = { ...
                sprintf('a1=%.6g, b1=%.6g', p.a1, p.b1), ...
                sprintf('a2=%.6g, b2=%.6g, c2=%.6g', p.a2, p.b2, p.c2), ...
                sprintf('psi_max=%.3g, X_max=%.3g, Ngrid=%d, shiftMin=%d', ...
                        view.psi_max, view.X_max, view.Ngrid, view.shiftMin) ...
            };
            drawnow limitrate;
        catch ME
            txtInfo.Value = { ...
                '[ERROR]', ME.message ...
            };
        end
    end

end

% ===================== local utility funcs =====================
function [sld, edt] = add_slider_row_(g, row, name, lims, val)
    uilabel(g,'Text',name,'HorizontalAlignment','right');
    sld = uislider(g,'Limits',lims,'Value',val);
    edt = uieditfield(g,'numeric','Value',val,'Limits',lims);
    sld.Layout.Row = row; sld.Layout.Column = 2;
    edt.Layout.Row = row; edt.Layout.Column = 3;
end

function set_slider_edit_(sld, edt, v)
    sld.Value = v;
    edt.Value = v;
end

function hook_slider_(sld, edt, setter)
    sld.ValueChangingFcn = @(src,ev)on_change(ev.Value);
    sld.ValueChangedFcn  = @(src,ev)on_change(src.Value);
    function on_change(v)
        edt.Value = v;
        setter(v);
    end
end

function hook_edit_(edt, sld, setter)
    edt.ValueChangedFcn = @(src,~)on_edit(src.Value);
    function on_edit(v)
        % clamp into slider limits
        v = min(max(v, sld.Limits(1)), sld.Limits(2));
        sld.Value = v;
        edt.Value = v;
        setter(v);
    end
end

function [psi, Fpsi_v, X, FX_v] = eval_curves_(p, view)
    % define polynomials
    Fpsi = @(psi) 0.5*p.a1*psi.^2 + (1/factorial(4))*p.b1*psi.^4;
    FX   = @(X)   0.5*p.a2*X.^2   + (1/factorial(4))*p.b2*X.^4 + (1/factorial(6))*p.c2*X.^6;

    N = max(101, round(view.Ngrid));
    psi = linspace(-view.psi_max, view.psi_max, N);
    X   = linspace(-view.X_max,   view.X_max,   N);

    Fpsi_v = Fpsi(psi);
    FX_v   = FX(X);

    if view.shiftMin
        Fpsi_v = Fpsi_v - min(Fpsi_v);
        FX_v   = FX_v   - min(FX_v);
    end
end

function s = view_ylabel_(shiftMin, which)
    if shiftMin
        if which=="psi"
            s = '$F_\psi(\psi)-\min$';
        else
            s = '$F_X(X)-\min$';
        end
    else
        if which=="psi"
            s = '$F_\psi(\psi)$';
        else
            s = '$F_X(X)$';
        end
    end
end

function t = title_psi_(p)
    t = sprintf('$F_\\psi=\\frac12 a_1\\psi^2+\\frac{1}{4!}b_1\\psi^4\\quad$');
end

function t = title_X_(p)
    t = sprintf('$F_X=\\frac12 a_2X^2+\\frac{1}{4!}b_2X^4+\\frac{1}{6!}c_2X^6\\quad$');
end

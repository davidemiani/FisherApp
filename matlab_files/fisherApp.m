classdef fisherApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                   matlab.ui.Figure
        Formula                    matlab.ui.control.Button
        TestButton                 matlab.ui.control.Button
        OasisPanel                 matlab.ui.container.Panel
        OasisKnobLabel             matlab.ui.control.Label
        OasisKnob                  matlab.ui.control.DiscreteKnob
        VEditFieldLabel1           matlab.ui.control.Label
        VEditField                 matlab.ui.control.NumericEditField
        VEditFieldLabel2           matlab.ui.control.Label
        OasisLamp                  matlab.ui.control.Lamp
        CoefficientsPanel          matlab.ui.container.Panel
        DSliderLabel2              matlab.ui.control.Label
        dSliderLabel               matlab.ui.control.Label
        DSlider                    matlab.ui.control.Slider
        ASliderLabel2              matlab.ui.control.Label
        aSliderLabel               matlab.ui.control.Label
        ASlider                    matlab.ui.control.Slider
        BSliderLabel2              matlab.ui.control.Label
        bSliderLabel               matlab.ui.control.Label
        BSlider                    matlab.ui.control.Slider
        GridPanel                  matlab.ui.container.Panel
        SpaceLengthEditFieldLabel  matlab.ui.control.Label
        SpaceEditField             matlab.ui.control.NumericEditField
        TimeLengthEditFieldLabel   matlab.ui.control.Label
        TimeEditField              matlab.ui.control.NumericEditField
        FocusedICSwitchLabel       matlab.ui.control.Label
        ICSwitch                   matlab.ui.control.Switch
        SpaceLabel2                matlab.ui.control.Label
        TimeLabel2                 matlab.ui.control.Label
        ResetValuesButton          matlab.ui.control.Button
    end


    properties (Access = public)
        Init
    end

    methods (Access = public)
        
        function oasisManager(app)
            switch app.OasisKnob.Value
                case 'Off'
                    app.OasisLamp.Color = [1,0,0];
                    app.VEditField.Enable = 'off';
                    app.VEditFieldLabel1.Enable = 'off';
                    app.VEditFieldLabel2.Enable = 'off';
                    app.DSliderLabel2.Text = 'e-1';
                    app.ASliderLabel2.Text = 'e-3';
                    app.BSliderLabel2.Text = 'e-6';
                otherwise
                    app.OasisLamp.Color = [0,1,0];
                    app.VEditField.Enable = 'on';
                    app.VEditFieldLabel1.Enable = 'on';
                    app.VEditFieldLabel2.Enable = 'on';
                    app.DSliderLabel2.Text = 'e-0';
                    app.ASliderLabel2.Text = 'e-2';
                    app.BSliderLabel2.Text = 'e-6';
            end
        end
        
        function L = fisherMat(~,n)
            L = -2*eye(n)+diag(ones(n-1,1),+1)+diag(ones(n-1,1),-1);
        end
        
        function uval = fisherWave(~,x,t,h,k,a,v,l)
            if l==0
                uval = a;
            else
                x = x(:);
                uval = zeros(size(x));
                iii = round(v*(t*k)/h+1);
                iif = round(iii+l/h);
                if iif>length(x)
                    iif = length(x);
                end
                hmi = iif-iii+1;
                uval(iii:iif) = a*ones(hmi,1);
            end
        end
        
    end


    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clc
            set(0,'DefaultFigurePosition',get(0,'screensize'))
            prop = properties(app);
            app.Init.Names = cell(0);
            app.Init.Values = cell(0);
            k = 1;
            for i=1:length(prop)
                if isprop(app.(prop{i}),'Value')
                    app.Init.Names{k} = prop{i};
                    app.Init.Values{k} = app.(prop{i}).Value;
                    k = k+1;
                end
            end
        end

        % Button pushed function: TestButton
        function TestButtonPushed(app, ~)
            app.UIFigure.Visible = 'off';
            
            x0 = 0;
            xf = app.SpaceEditField.Value;
            m = xf*10+1;
            x = linspace(x0,xf,m)';
            h = x(2)-x(1);
            
            t0 = 0;
            tf = app.TimeEditField.Value;
            n = tf*10+1;
            t = linspace(t0,tf,n);
            k = t(2)-t(1);
            
            c = zeros(m,n);
            
            if strcmp(app.OasisKnob.Value,'Off')
                d = app.DSlider.Value*10^-1;
                a = app.ASlider.Value*10^-3;
                b = app.BSlider.Value*10^-6;
                v = 0;
                l = 0;
            else
                d = app.DSlider.Value;
                a = app.ASlider.Value*10^-2;
                b = app.BSlider.Value*10^-5;
                v = app.VEditField.Value;
                switch app.OasisKnob.Value
                    case 'Low'
                        l = 1;
                    case 'Medium'
                        l = 3;
                    case 'High'
                        l = 5;
                end
            end
            r = d*k/(h^2);
            
            if strcmp(app.ICSwitch.Value,'Off')
                ss0 = x0;
                ssf = xf;
            else
                ss0 = x0;
                ssf = 5;
            end
            hms = ssf-ss0;
            ii0 = round(ss0/h+1);
            iif = round(ssf/h+1);
            X = x(ii0:iif);
            c(ii0:iif,1) = 5000*exp(0.1*X).*sin(pi*X/hms);
            
            
            A = eye(m-2)+r*fisherMat(app,m-2);
            B = (eye(m-2)-r*fisherMat(app,m-2))^-1;
            for j=1:n-1
               C = c(2:m-1,j);
               X = x(2:m-1,1);
               c(2:m-1,j+1) = B*A*C+B*(fisherWave(app,X,j,h,k,a,v,l).*C)-b*B*C.^2;
            end
            
            h = figure;
            hold on
            colormap jet
            shading flat
            grid on
            surf(t,x,c,'EdgeColor','none','FaceLighting','gouraud')
            for j=1:round(n/(n/5)):n
               plot3(t(j)*ones(size(x)),x,c(:,j),'k')
            end
            for i=1:round(m/(m/5)):m
                plot3(t,x(i)*ones(size(t)),c(i,:),'k')
            end
            view(3)
            set(gca,'FontSize',30)
            xlabel('time [min]')
            ylabel('space [cm]')
            zlabel('bacteria conc. [cm^{-1}]')
            nb0 = round(trapz(x,c(:,1))*10^-3);
            nbf = round(trapz(x,c(:,end))*10^-3);
            dfb = (nbf-nb0)/nb0*100;
            title({['d=',num2str(d),'   ', ...
                    'a=',num2str(a),'   ', ...
                    'b=',num2str(b),'   ', ...
                    'v=',num2str(v)]     ; ...
                   ['n_b_0=',num2str(nb0),'k   ', ...
                    'n_b_f=',num2str(nbf),'k   ', ...
                    '\Delta_%=',num2str(dfb),'%']})
            uiwait(h)
            
            app.UIFigure.Visible = 'on';
        end

        % Value changed function: OasisKnob
        function OasisKnobValueChanged(app, ~)
            oasisManager(app)
        end

        % Button pushed function: ResetValuesButton
        function ResetValuesButtonPushed(app, ~)
            names = app.Init.Names;
            values = app.Init.Values;
            for i=1:length(names)
                app.(names{i}).Value = values{i};
            end
            oasisManager(app)
        end

        % Button pushed function: Formula
        function FormulaButtonPushed(app, ~)
            line1 = ['Fisher proposed this equation in 1937 to describe the spatial ', ...
                     'spread of an advantageous allele and explore its wave solution', ...
                     's. In this case, we can use it to describe bacterial culture e', ...
                     'volution.'];
            line2 = ['With this application, you can easly set the equation coeffici', ...
                     'ents and see how solution change:'];
            line3 = '   - d: diffusion parameter';
            line4 = '   - a: growth parameter';
            line5 = '   - b: resource constraints parameter';
            line6 = '   - v: oasis speed';
            line7 = 'You can find more information on attached paper. Good testing!';
            msg = {line1;blanks(1);line2;line3;line4;line5;line6;blanks(1);line7};
            uialert(app.UIFigure,msg,'Equation''s Info','icon','info');
        end

        % Value changed function: DSlider
        function DSliderValueChanged(app, ~)
            app.DSlider.Value = round(app.DSlider.Value,1);
        end

        % Value changed function: ASlider
        function ASliderValueChanged(app, ~)
            app.ASlider.Value = round(app.ASlider.Value,1);
        end

        % Value changed function: BSlider
        function BSliderValueChanged(app, ~)
            app.BSlider.Value = round(app.BSlider.Value,1);
        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 649 385];
            app.UIFigure.Name = '1D Bacterial Culture Simulation';
            app.UIFigure.Resize = 'off';
            setAutoResize(app, app.UIFigure, true)

            % Create Formula
            app.Formula = uibutton(app.UIFigure, 'push');
            app.Formula.ButtonPushedFcn = createCallbackFcn(app, @FormulaButtonPushed, true);
            app.Formula.Icon = 'fisherEq.png';
            app.Formula.BackgroundColor = [0.9373 0.9373 0.9373];
            app.Formula.Position = [23 320 605 55];
            app.Formula.Text = '';

            % Create TestButton
            app.TestButton = uibutton(app.UIFigure, 'push');
            app.TestButton.ButtonPushedFcn = createCallbackFcn(app, @TestButtonPushed, true);
            app.TestButton.Icon = 'matlab_24.png';
            app.TestButton.IconAlignment = 'top';
            app.TestButton.Position = [408 34 88 55];
            app.TestButton.Text = 'Test';

            % Create OasisPanel
            app.OasisPanel = uipanel(app.UIFigure);
            app.OasisPanel.TitlePosition = 'centertop';
            app.OasisPanel.Title = 'OASIS SIMULATION';
            app.OasisPanel.Position = [397 122 209 187];

            % Create OasisKnobLabel
            app.OasisKnobLabel = uilabel(app.OasisPanel);
            app.OasisKnobLabel.HorizontalAlignment = 'center';
            app.OasisKnobLabel.Position = [37 58 94 15];
            app.OasisKnobLabel.Text = 'Oasis Thickness';

            % Create OasisKnob
            app.OasisKnob = uiknob(app.OasisPanel, 'discrete');
            app.OasisKnob.ValueChangedFcn = createCallbackFcn(app, @OasisKnobValueChanged, true);
            app.OasisKnob.Position = [72 80 60 60];

            % Create VEditFieldLabel1
            app.VEditFieldLabel1 = uilabel(app.OasisPanel);
            app.VEditFieldLabel1.Enable = 'off';
            app.VEditFieldLabel1.HorizontalAlignment = 'right';
            app.VEditFieldLabel1.Position = [34 15 22 15];
            app.VEditFieldLabel1.Text = 'v';

            % Create VEditField
            app.VEditField = uieditfield(app.OasisPanel, 'numeric');
            app.VEditField.Limits = [0 Inf];
            app.VEditField.Enable = 'off';
            app.VEditField.Position = [64 12 77 22];
            app.VEditField.Value = 1;

            % Create VEditFieldLabel2
            app.VEditFieldLabel2 = uilabel(app.OasisPanel);
            app.VEditFieldLabel2.Enable = 'off';
            app.VEditFieldLabel2.Position = [146 15 46 15];
            app.VEditFieldLabel2.Text = 'cm/min';

            % Create OasisLamp
            app.OasisLamp = uilamp(app.OasisPanel);
            app.OasisLamp.Position = [150 55 20 20];
            app.OasisLamp.Color = [1 0 0];

            % Create CoefficientsPanel
            app.CoefficientsPanel = uipanel(app.UIFigure);
            app.CoefficientsPanel.TitlePosition = 'centertop';
            app.CoefficientsPanel.Title = 'COEFFICIENTS';
            app.CoefficientsPanel.Position = [45 122 330 187];

            % Create DSliderLabel2
            app.DSliderLabel2 = uilabel(app.CoefficientsPanel);
            app.DSliderLabel2.Position = [295 141 25 15];
            app.DSliderLabel2.Text = 'e-1';

            % Create dSliderLabel
            app.dSliderLabel = uilabel(app.CoefficientsPanel);
            app.dSliderLabel.HorizontalAlignment = 'right';
            app.dSliderLabel.Position = [-4 141 25 15];
            app.dSliderLabel.Text = 'd';

            % Create DSlider
            app.DSlider = uislider(app.CoefficientsPanel);
            app.DSlider.Limits = [0 5];
            app.DSlider.MajorTicks = [0 1 2 3 4 5];
            app.DSlider.ValueChangedFcn = createCallbackFcn(app, @DSliderValueChanged, true);
            app.DSlider.Position = [42 147 238 3];
            app.DSlider.Value = 1;

            % Create ASliderLabel2
            app.ASliderLabel2 = uilabel(app.CoefficientsPanel);
            app.ASliderLabel2.Position = [295 91 25 15];
            app.ASliderLabel2.Text = 'e-3';

            % Create aSliderLabel
            app.aSliderLabel = uilabel(app.CoefficientsPanel);
            app.aSliderLabel.HorizontalAlignment = 'right';
            app.aSliderLabel.Position = [-4 91 25 15];
            app.aSliderLabel.Text = 'a';

            % Create ASlider
            app.ASlider = uislider(app.CoefficientsPanel);
            app.ASlider.Limits = [0 5];
            app.ASlider.MajorTicks = [0 1 2 3 4 5];
            app.ASlider.ValueChangedFcn = createCallbackFcn(app, @ASliderValueChanged, true);
            app.ASlider.Position = [42 97 238 3];

            % Create BSliderLabel2
            app.BSliderLabel2 = uilabel(app.CoefficientsPanel);
            app.BSliderLabel2.Position = [295 41 25 15];
            app.BSliderLabel2.Text = 'e-6';

            % Create bSliderLabel
            app.bSliderLabel = uilabel(app.CoefficientsPanel);
            app.bSliderLabel.HorizontalAlignment = 'right';
            app.bSliderLabel.Position = [-4 41 25 15];
            app.bSliderLabel.Text = 'b';

            % Create BSlider
            app.BSlider = uislider(app.CoefficientsPanel);
            app.BSlider.Limits = [0 5];
            app.BSlider.MajorTicks = [0 1 2 3 4 5];
            app.BSlider.ValueChangedFcn = createCallbackFcn(app, @BSliderValueChanged, true);
            app.BSlider.Position = [42 47 238 3];

            % Create GridPanel
            app.GridPanel = uipanel(app.UIFigure);
            app.GridPanel.TitlePosition = 'centertop';
            app.GridPanel.Title = 'GRID DEFINITION';
            app.GridPanel.Position = [45 12 330 99];

            % Create SpaceLengthEditFieldLabel
            app.SpaceLengthEditFieldLabel = uilabel(app.GridPanel);
            app.SpaceLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.SpaceLengthEditFieldLabel.Position = [23.03125 50 81 15];
            app.SpaceLengthEditFieldLabel.Text = 'Space Length';

            % Create SpaceEditField
            app.SpaceEditField = uieditfield(app.GridPanel, 'numeric');
            app.SpaceEditField.Limits = [20 30];
            app.SpaceEditField.RoundFractionalValues = 'on';
            app.SpaceEditField.ValueDisplayFormat = '%.0f';
            app.SpaceEditField.Position = [111.03125 46 54.953125 22];
            app.SpaceEditField.Value = 20;

            % Create TimeLengthEditFieldLabel
            app.TimeLengthEditFieldLabel = uilabel(app.GridPanel);
            app.TimeLengthEditFieldLabel.HorizontalAlignment = 'right';
            app.TimeLengthEditFieldLabel.Position = [31.03125 17 73 15];
            app.TimeLengthEditFieldLabel.Text = 'Time Length';

            % Create TimeEditField
            app.TimeEditField = uieditfield(app.GridPanel, 'numeric');
            app.TimeEditField.Limits = [10 240];
            app.TimeEditField.RoundFractionalValues = 'on';
            app.TimeEditField.Position = [111.03125 13 54.953125 22];
            app.TimeEditField.Value = 20;

            % Create FocusedICSwitchLabel
            app.FocusedICSwitchLabel = uilabel(app.GridPanel);
            app.FocusedICSwitchLabel.HorizontalAlignment = 'center';
            app.FocusedICSwitchLabel.Position = [229 15 67 15];
            app.FocusedICSwitchLabel.Text = 'Focused IC';

            % Create ICSwitch
            app.ICSwitch = uiswitch(app.GridPanel, 'slider');
            app.ICSwitch.Position = [240 45 45 20];

            % Create SpaceLabel2
            app.SpaceLabel2 = uilabel(app.GridPanel);
            app.SpaceLabel2.Position = [173 50 25 15];
            app.SpaceLabel2.Text = 'cm';

            % Create TimeLabel2
            app.TimeLabel2 = uilabel(app.GridPanel);
            app.TimeLabel2.Position = [173 17 25 15];
            app.TimeLabel2.Text = 'min';

            % Create ResetValuesButton
            app.ResetValuesButton = uibutton(app.UIFigure, 'push');
            app.ResetValuesButton.ButtonPushedFcn = createCallbackFcn(app, @ResetValuesButtonPushed, true);
            app.ResetValuesButton.Icon = 'cancel_24.png';
            app.ResetValuesButton.IconAlignment = 'top';
            app.ResetValuesButton.Position = [507 34 90.5 55];
            app.ResetValuesButton.Text = 'Reset Values';
        end
    end

    methods (Access = public)

        % Construct app
        function app = fisherApp()

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
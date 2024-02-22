% Yield calculation
annewmodel=iAQY970
s=optimizeCbModel(nnewmodel,'max','one');

% Find the metabolite of ATP
nad_index = find(strcmp(annewmodel.mets, 'atp[c]'));

% Find reactions involved atp[c]
atp_production_reaction = find(annewmodel.S(nad_index, :) > 0);

% Obtain the flux values of these reactions
atp_production_flux = s.v(atp_production_reaction);

% Calculate ATP yield
atp_yield = sum(atp_production_flux);

fprintf('ATP Production Rate: %.4f mmol/gDW/h\n', atp_yield);

% Find the metabolite of NADPH
nadph_index = find(strcmp(annewmodel.mets, 'nadph[c]'));

% Find reactions involved nadph[c]
nadph_production_reaction = find(annewmodel.S(nadph_index, :) > 0);

% Obtain the flux values of these reactions
nadph_production_flux = s.v(nadph_production_reaction);

% Calculate NADPH yield
nadph_yield = sum(nadph_production_flux);

fprintf('NADPH Production Rate: %.4f mmol/gDW/h\n', nadph_yield);

% Find the metabolite of NADH
nadh_index = find(strcmp(annewmodel.mets, 'nadh[c]'));

% Find reactions involved nadh[c]
nadh_production_reaction = find(annewmodel.S(nadh_index, :) > 0);

% Obtain the flux values of these reactions
nadh_production_flux = s.v(nadh_production_reaction);

% Calculate NADH yield
nadh_yield = sum(nadh_production_flux);

fprintf('NADH Production Rate: %.4f mmol/gDW/h\n', nadh_yield);

% Find the metabolite of fixed NH3
fixedNH3_index = find(strcmp(annewmodel.mets, 'fixedNH3[c]'));

% Find reactions involved fixedNH3[c]
fixedNH3_production_reaction = find(annewmodel.S(fixedNH3_index, :) > 0);

% Obtain the flux values of these reactions
fixedNH3_production_flux = s.v(fixedNH3_production_reaction);

% Calculate fixed NH3 yield
fixedNH3_yield = sum(fixedNH3_production_flux);

fprintf('fixedNH3 Production Rate: %.4f mmol/gDW/h\n', fixedNH3_yield);


function d = entry_margin_diagnostics()
global ce cw fspacee fspacew smine smaxe sminw smaxw INV_ME INV_MW INV_AG ...
    egrid ygrid k Ny Perg piy theta phi gamma alpha eta W Rnew qu pu F ...
    lambda_u productivity_weight r delta kappa zeta cxw

d = struct();

ag = INV_AG(:);
Na = numel(ag);
Ns = Na*k*Ny;
if numel(INV_ME) ~= Ns
    error('entry_margin_diagnostics: INV_ME size %d ~= Na*k*Ny = %d', numel(INV_ME), Ns);
end

S   = gridmake(ag, (1:k)', (1:Ny)');
a_s = S(:,1);
z_s = S(:,2);
y_s = S(:,3);

me   = INV_ME(:);
mw   = INV_MW(:);
mtot = me + mw;
Mass = sum(mtot);
me   = me/Mass;
mw   = mw/Mass;
mtot = mtot/Mass;

aw_eval  = min(max(a_s, sminw(1)), smaxw(1));
ae_eval  = min(max(a_s, smine(1)), smaxe(1));
VW       = funeval_fast(cw(:,1), fspacew, aw_eval, z_s, y_s);
VE       = funeval_fast(ce(:,1), fspacee, ae_eval, z_s, y_s);
prefer_E = (VE >= VW);

astar_zy = nan(k, Ny);
for iz = 1:k
    for iy = 1:Ny
        sel = (z_s == iz) & (y_s == iy);
        pe  = prefer_E(sel);
        aa  = a_s(sel);
        j   = find(pe, 1, 'first');
        if isempty(j)
            astar_zy(iz,iy) = Inf;
        else
            astar_zy(iz,iy) = aa(j);
        end
    end
end

piy_c   = piy(:);
piy_c   = piy_c/sum(piy_c);
astar_z = nan(k,1);
for iz = 1:k
    v   = astar_zy(iz,:)';
    fin = isfinite(v);
    if any(fin)
        wgt         = piy_c(fin)/sum(piy_c(fin));
        astar_z(iz) = sum(v(fin).*wgt);
    end
end
d.astar_zy = astar_zy;
d.astar_z  = astar_z;
d.zP       = exp(egrid(:));
d.pi_z     = Perg(:);
d.frac_z_never_enter = sum(Perg(:) .* all(~isfinite(astar_zy),2));

[a_sorted, ord] = sort(a_s);
cw_mass       = cumsum(mtot(ord));
median_wealth = a_sorted(find(cw_mass >= 0.50, 1, 'first'));
p25_wealth    = a_sorted(find(cw_mass >= 0.25, 1, 'first'));
p75_wealth    = a_sorted(find(cw_mass >= 0.75, 1, 'first'));
d.median_wealth = median_wealth;
d.p25_wealth    = p25_wealth;
d.p75_wealth    = p75_wealth;
d.mean_wealth   = sum(mtot .* a_s);

cz            = cumsum(Perg(:));
top10_z_nodes = find(cz > 0.90);
top1_z_nodes  = find(cz > 0.99);
is_top10_z    = ismember(z_s, top10_z_nodes);
is_top1_z     = ismember(z_s, top1_z_nodes);
is_poor       = (a_s <= median_wealth);

safe = @(num,den) num/max(den,1e-14);
d.ent_share_overall        = sum(me);
d.ent_top_abil_decile      = safe(sum(me(is_top10_z)),        sum(mtot(is_top10_z)));
d.ent_top_abil_pct1        = safe(sum(me(is_top1_z)),         sum(mtot(is_top1_z)));
d.ent_bottom_abil_90       = safe(sum(me(~is_top10_z)),       sum(mtot(~is_top10_z)));
d.poor_talented_ent_rate   = safe(sum(me(is_poor & is_top10_z)),  sum(mtot(is_poor & is_top10_z)));
d.rich_talented_ent_rate   = safe(sum(me(~is_poor & is_top10_z)), sum(mtot(~is_poor & is_top10_z)));
d.poor_talented_mass       = sum(mtot(is_poor & is_top10_z));
d.ent_share_bottom_half_wealth = safe(sum(me(is_poor)), sum(me));
d.ent_mass_poor_talented   = sum(me(is_poor & is_top10_z));

zlev   = exp(egrid(:));
z_of_s = zlev(z_s);
d.talent_utilisation = safe(sum(me .* z_of_s), sum(mtot .* z_of_s));

astar_top = astar_z(top10_z_nodes);
d.astar_top_abil = median(astar_top(isfinite(astar_top)));
d.astar_top_abil_pctile = safe(sum(mtot(a_s <= d.astar_top_abil)), 1);

Wealth_E = sum(me .* a_s);
Wealth_W = sum(mw .* a_s);
d.wealth_ent = Wealth_E;
d.wealth_wrk = Wealth_W;
d.ent_wealth_share = safe(Wealth_E, Wealth_E + Wealth_W);
d.one_over_theta   = 1/theta;
d.wealth_share_gap = d.ent_wealth_share - 1/theta;
d.mean_wealth_ent  = safe(Wealth_E, sum(me));
d.mean_wealth_wrk  = safe(Wealth_W, sum(mw));

amin_      = min(ag);
near_floor = (a_s <= amin_ + 1e-6 + 0.002*(max(ag)-amin_));
d.worker_mass_at_floor = safe(sum(mw(near_floor)), sum(mw));
d.ent_mass_at_floor    = safe(sum(me(near_floor)), sum(me));
d.mass_at_amax         = safe(sum(mtot(a_s >= max(ag)*0.999)), 1);

sel_e = me > 1e-14;
if any(sel_e)
    A_e = a_s(sel_e);
    Z_e = productivity_weight .* exp(egrid(z_s(sel_e)));
    m_e = me(sel_e);
    [Kn_e, Ku_e, Kagg_e, L_e, Y_e, D_e, mu_e, bind_e] = compute_capital_allocation( ...
        A_e, Z_e, Rnew, qu, phi, gamma, alpha, eta, W, F, pu, theta, 'vectorized', lambda_u);
    wgt = m_e/sum(m_e);
    d.binding_share = sum(wgt .* double(bind_e));
    d.mean_mu       = sum(wgt .* mu_e);

    val_n = Kn_e;
    val_u = pu.*Ku_e;
    d.used_value_share = safe(sum(m_e.*val_u), sum(m_e.*(val_n+val_u)));

    poor_e = A_e <= median_wealth;
    d.used_value_share_poor = safe(sum(m_e(poor_e).*val_u(poor_e)), sum(m_e(poor_e).*(val_n(poor_e)+val_u(poor_e))));
    d.used_value_share_rich = safe(sum(m_e(~poor_e).*val_u(~poor_e)), sum(m_e(~poor_e).*(val_n(~poor_e)+val_u(~poor_e))));
    d.binding_share_poor = safe(sum(m_e(poor_e).*double(bind_e(poor_e))), sum(m_e(poor_e)));
    d.binding_share_rich = safe(sum(m_e(~poor_e).*double(bind_e(~poor_e))), sum(m_e(~poor_e)));

    mpk  = (1-alpha)*eta .* Y_e ./ max(Kagg_e,1e-12);
    good = isfinite(mpk) & mpk > 0;
    lw   = m_e(good)/sum(m_e(good));
    lm   = log(mpk(good));
    d.mpk_mean_log = sum(lw.*lm);
    d.mpk_sd_log   = sqrt(max(sum(lw.*(lm - d.mpk_mean_log).^2), 0));

    Yagg     = sum(m_e.*Y_e);
    Lagg     = sum(m_e.*L_e);
    Kagg_tot = sum(m_e.*Kagg_e);
    d.Y_agg = Yagg;
    d.L_agg = Lagg;
    d.K_agg = Kagg_tot;
    d.TFP_agg = safe(Yagg, (Lagg^alpha * Kagg_tot^(1-alpha))^eta);
    d.n_firms = sum(m_e);

    d.avg_ability_ent = safe(sum(m_e.*Z_e), sum(m_e));
    d.avg_ability_all = safe(sum(mtot.*exp(egrid(z_s))), 1);
end

entrants_from_W = sum(mw .* double(prefer_E & (a_s > smine(1))));
d.entrants_from_wagework_mass = entrants_from_W;
d.entry_rate_from_wagework    = safe(entrants_from_W, sum(mw));

iy_med = ceil(Ny/2);
a0     = min(ag);
maxT   = 400;
yte    = nan(k,1);
for iz = 1:k
    target = astar_zy(iz, iy_med);
    if ~isfinite(target), continue; end
    aa = a0;
    if aa >= target, yte(iz) = 0; continue; end
    for t = 1:maxT
        ae    = min(max(aa, sminw(1)), smaxw(1));
        anext = funeval_fast(cxw, fspacew, ae, iz, iy_med);
        anext = min(max(anext, sminw(1)), smaxw(1));
        if anext <= aa + 1e-12
            yte(iz) = Inf;
            break
        end
        aa = anext;
        if aa >= target, yte(iz) = t; break; end
    end
    if isnan(yte(iz)), yte(iz) = Inf; end
end
d.years_to_entry_z = yte;
yt_top = yte(top10_z_nodes);
d.years_to_entry_top_abil    = median(yt_top);
d.frac_top_abil_ladder_broken = mean(~isfinite(yt_top));
finite_top = yt_top(isfinite(yt_top));
if ~isempty(finite_top)
    d.years_to_entry_top_abil_finite = median(finite_top);
else
    d.years_to_entry_top_abil_finite = Inf;
end

a_floor_next = funeval_fast(cxw, fspacew, min(max(a0,sminw(1)),smaxw(1)), ...
    max(top10_z_nodes(1),1), iy_med);
d.floor_worker_savings = a_floor_next - a0;
d.beta_gross_return    = NaN;

fprintf('  [entry diag] years-to-entry (top abil, median y) = %.1f | ladder broken for %.0f%% of top nodes\n', ...
    d.years_to_entry_top_abil, 100*d.frac_top_abil_ladder_broken);
fprintf('  [entry diag] top-abil-decile ent=%.4f | poor&talented ent=%.4f | ent wealth share=%.4f (1/theta=%.4f)\n', ...
    d.ent_top_abil_decile, d.poor_talented_ent_rate, d.ent_wealth_share, 1/theta);
fprintf('  [entry diag] a*(top abil)=%.4f at wealth pctile %.3f | median wealth=%.4f | workers at floor=%.3f\n', ...
    d.astar_top_abil, d.astar_top_abil_pctile, d.median_wealth, d.worker_mass_at_floor);
end

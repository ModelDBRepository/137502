function checkStriatumInputs(SIMPARAMS)

if ~isfloat(SIMPARAMS.sim.tstart); error('tstart is not a float'); end
if ~isfloat(SIMPARAMS.sim.tfinal); error('tfinal is not a float'); end
if ~isfloat(SIMPARAMS.sim.dt); error('dt is not a float'); end
if ~isfloat(SIMPARAMS.physiology.MSparams); error('MSparams is not a float'); end
if ~isfloat(SIMPARAMS.physiology.FSparams); error('FSparams is not a float'); end
if ~isfloat(SIMPARAMS.physiology.Eglu); error('Eglu is not a float'); end
if ~isfloat(SIMPARAMS.physiology.Egaba); error('Egaba is not a float'); end
if ~isfloat(SIMPARAMS.physiology.ts_glu_AMPA); error('ts_glu_AMPA is not a float'); end
if ~isfloat(SIMPARAMS.physiology.ts_glu_NMDA); error('ts_glu_NMDA is not a float'); end
if ~isfloat(SIMPARAMS.physiology.ts_gaba); error('ts_gaba is not a float'); end
if ~isfloat(SIMPARAMS.physiology.tau_fsgap); error('tau_fsgap is not a float'); end
if ~isinteger(SIMPARAMS.sim.MSspikebuffer); error('MSspikebuffer is not a integer'); end
if ~isinteger(SIMPARAMS.sim.FSspikebuffer); error('FSspikebuffer is not a integer'); end
if ~isfloat(SIMPARAMS.sim.initVms); error('initVms is not a float'); end
if ~isfloat(SIMPARAMS.sim.initUms); error('initUms is not a float'); end
if ~isfloat(SIMPARAMS.sim.initVfs); error('initVfs is not a float'); end
if ~isfloat(SIMPARAMS.sim.initUfs); error('initUfs is not a float'); end
if ~isfloat(SIMPARAMS.sim.initVgapfs); error('initVgapfs is not a float'); end
if ~isfloat(SIMPARAMS.sim.SpikeEventQue_MS); error('SpikeEventQue_MS is not a float'); end
if ~isfloat(SIMPARAMS.sim.SpikeEventQue_FS); error('SpikeEventQue_FS is not a float'); end
if ~isfloat(SIMPARAMS.initCTX); error('initCTX is not a float'); end
if ~isfloat(SIMPARAMS.sim.Iinj_MS); error('Iinj is not a float'); end
if ~isfloat(SIMPARAMS.sim.Iinj_FS); error('Iinj is not a float'); end
if ~isinteger(SIMPARAMS.net.Cctms); error('Cctms is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cctms_b); error('Cctms_b is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cctms_d); error('Cctms_d is not a integer'); end
if ~isfloat(SIMPARAMS.net.Cctms_w); error('Cctms_w is not a float'); end
if ~isinteger(SIMPARAMS.net.Cmsms); error('Cmsms is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cmsms_b); error('Cmsms_b is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cmsms_d); error('Cmsms_d is not a integer'); end
if ~isfloat(SIMPARAMS.net.Cmsms_w); error('Cmsms_w is not a float'); end
if ~isinteger(SIMPARAMS.net.Cfsms); error('Cfsms is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cfsms_b); error('Cfsms_b is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cfsms_d); error('Cfsms_d is not a integer'); end
if ~isfloat(SIMPARAMS.net.Cfsms_w); error('Cfsms_w is not a float'); end
if ~isinteger(SIMPARAMS.net.Cctfs); error('Cctfs is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cctfs_b); error('Cctfs_b is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cctfs_d); error('Cctfs_d is not a integer'); end
if ~isfloat(SIMPARAMS.net.Cctfs_w); error('Cctfs_w is not a float'); end
if ~isinteger(SIMPARAMS.net.Cfsfs); error('Cfsfs is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cfsfs_b); error('Cfsfs_b is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cfsfs_d); error('Cfsfs_d is not a integer'); end
if ~isfloat(SIMPARAMS.net.Cfsfs_w); error('Cfsfs_w is not a float'); end
if ~isinteger(SIMPARAMS.net.Cgapfs); error('Cgapfs is not a integer'); end
if ~isinteger(SIMPARAMS.net.Cgapfs_b); error('Cgapfs_b is not a integer'); end
if ~isfloat(SIMPARAMS.net.Cgapfs_w); error('Cgapfs_w is not a float'); end
if ~isinteger(SIMPARAMS.net.Pgapfs); error('Pgapfs is not a integer'); end
if ~isfloat(SIMPARAMS.CTX_state); error('CTX_state is not a float'); end
if ~isinteger(SIMPARAMS.net.CHAN1_MS); error('CHAN1_MS is not a integer'); end
if ~isinteger(SIMPARAMS.net.CHAN1_FS); error('CHAN1_FS is not a integer'); end
if ~isinteger(SIMPARAMS.net.CHAN2_MS); error('CHAN2_MS is not a integer'); end
if ~isinteger(SIMPARAMS.net.CHAN2_FS); error('CHAN2_FS is not a integer'); end
if ~isinteger(SIMPARAMS.input.CTX.N_MSSEG); error('N_MSSEG is not a integer'); end
if ~isfloat(SIMPARAMS.input.CTX.r_MSSEG); error('r_MSSEG is not a float'); end
if ~isfloat(SIMPARAMS.input.CTX.alpha_MSSEG); error('alpha_MSSEG is not a float'); end
if ~isinteger(SIMPARAMS.input.CTX.N_FSSEG); error('N_FSSEG is not a integer'); end
if ~isfloat(SIMPARAMS.input.CTX.r_FSSEG); error('r_FSSEG is not a float'); end
if ~isfloat(SIMPARAMS.input.CTX.alpha_FSSEG); error('alpha_FSSEG is not a float'); end
if ~isfloat(SIMPARAMS.physiology.glu_ratio); error('glu_ratio is not a float'); end
if ~isfloat(SIMPARAMS.physiology.DA); error('DA is not a float'); end
if ~isinteger(SIMPARAMS.sim.RecordChan_MS); error('RecordChan_MS is not a integer'); end
if ~isfloat(SIMPARAMS.input.PULSE.P); error('PULSE.P is not a float'); end
if ~isfloat(SIMPARAMS.input.PULSE.Nctx_ms); error('PULSE.Nctx_ms is not a float'); end
if ~isfloat(SIMPARAMS.input.PULSE.Nctx_fs); error('PULSE.Nctx_fs is not a float'); end
if ~isfloat(SIMPARAMS.input.PULSE.ts_spks); error('PULSE.ts_spks is not a float'); end
if ~isinteger(SIMPARAMS.sim.RANDSEED); error('RANDSEED is not a integer'); end



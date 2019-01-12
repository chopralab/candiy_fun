# Find tryptamines
obgrep 'c1ccc2c(c1)c(c[nH]2)CCN' all.smi > tryp.smi
obgrep -v 'c1ccc2c(c1)c(c[nH]2)CCN' all.smi > not_tryp.smi

# Find cathinones
obgrep 'O=C(c1ccccc1)C(N)C' not_tryp.smi > cath.smi
obgrep -v 'O=C(c1ccccc1)C(N)C' not_tryp.smi > not_tryp_cath.smi

# Find amphetamines
obgrep 'c1ccccc1[CH2]C(N)[CH3]' not_tryp_cath.smi > amph.smi
obgrep -v 'c1ccccc1[CH2]C(N)[CH3]' not_tryp_cath.smi > not_tryp_cath_amph.smi

# Find remaining phenethylamines
obgrep 'c1ccccc1CCN' not_tryp_cath_amph.smi > phen.smi
obgrep -v 'c1ccccc1CCN' not_tryp_cath_amph.smi > not_tryp_cath_amph_phen.smi

# Find cannabinols
obgrep 'CCCCCc1cc(O)c(C)cc1' not_tryp_cath_amph.smi > canb.smi
obgrep -v 'CCCCCc1cc(O)ccc1' not_tryp_cath_amph_phen.smi > other.smi


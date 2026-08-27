# axion-field object source files
src += os-axion-field-diag-define.f03 os-axion-field-define.f03 \
       os-axion-field-boundary.f03 os-axion-field-core.f03 \
       os-axion-field-solver.f03 os-axion-field-diagnostics.f03 \
       os-axion-field.f03 os-simulation-axion-hybrid.f03

# Explicit ordering for the split definition/implementation units
os-axion-field-define.f03: os-axion-field-diag-define.f03
os-axion-field-boundary.f03: os-axion-field-define.f03
os-axion-field-core.f03: os-axion-field-define.f03 os-axion-field-boundary.f03
os-axion-field-solver.f03: os-axion-field-define.f03 os-axion-field-boundary.f03
os-axion-field-diagnostics.f03: os-axion-field-diag-define.f03 os-axion-field-define.f03 os-axion-field-solver.f03

# Add local path to vpath
vpath %.f03 axion_field

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	@$(MAKE) -C ..

clean:
	@$(MAKE) -C .. clean

endif

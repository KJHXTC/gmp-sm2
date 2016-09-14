RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include src/subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif


# All Target
all: yl_sm2_gmp

# Tool invocations
yl_sm2_gmp: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc -D_DEBUG -L"../lib" -o"yl_sm2_gmp" $(OBJS) $(USER_OBJS) $(LIBS) -I ../include 
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) yl_sm2_gmp
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


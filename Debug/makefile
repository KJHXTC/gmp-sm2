RM := rm -rf
CC := gcc

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
all : demo liba

# Tool invocations
demo: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	$(CC) -o "demo" $(OBJS) $(USER_OBJS) $(LIBS) -I ../include -L"../lib" 
	@echo 'Finished building target: $@'
	@echo ' '

liba: 
	-$(AR) x ../lib/libgmp.a 
	-$(AR) -r libsm2enc.a $(OBJS) *.o
	-$(RM) *.o

# Other Targets
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) demo libsm2enc.a
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:


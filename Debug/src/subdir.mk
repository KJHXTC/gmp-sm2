################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/sm3.c \
../src/base_tools.c \
../src/ec_operations.c \
../src/int_arithmetic.c \
../src/poly_arithmetic.c \
../src/sm2_func.c \
../src/demo.c 

OBJS += \
./src/sm3.o \
./src/base_tools.o \
./src/ec_operations.o \
./src/int_arithmetic.o \
./src/poly_arithmetic.o \
./src/sm2_func.o \
./src/demo.o 

C_DEPS += \
./src/sm3.d \
./src/base_tools.d \
./src/ec_operations.d \
./src/int_arithmetic.d \
./src/poly_arithmetic.d \
./src/sm2_func.d \
./src/demo.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -I"../include" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



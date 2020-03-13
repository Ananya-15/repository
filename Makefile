.PHONY: clean All

All:
	@echo "----------Building project:[ HPCcoursework - Debug ]----------"
	@cd "HPCcoursework" && "$(MAKE)" -f  "HPCcoursework.mk"
clean:
	@echo "----------Cleaning project:[ HPCcoursework - Debug ]----------"
	@cd "HPCcoursework" && "$(MAKE)" -f  "HPCcoursework.mk" clean

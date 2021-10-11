import sys

def hello_world():
	print("Hello world!")

def good_morning():
	print("Good morning world!")

def good_night():
	print("Good night world!")

def main():
	function_names = sys.argv
	functions = {'hello_world': hello_world,
				'good_morning': good_morning,
				'good_night': good_night}
	if "all" in function_names:
		print("Executing entire workflow:")
		for function_name in functions:
			functions[function_name]()
	else:
		print(f"Executing functions {', '.join(sys.argv[1:])}:")
		for function_name in sys.argv[1:]:
			functions[function_name]()

main()

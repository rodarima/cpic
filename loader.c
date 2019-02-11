#include <stdio.h>
#include <assert.h>
#include <pthread.h>
#include <nanos6.h>
#include <nanos6/library-mode.h>
#include <nanos6/bootstrap.h>
#include <nanos6/api-check.h>


nanos6_api_versions_t const __user_code_expected_nanos6_api_versions = {
	.api_check_api_version = nanos6_api_check_api,
	.major_api_version = nanos6_major_api,
	
	.blocking_api_version = nanos6_blocking_api,
	.bootstrap_api_version = nanos6_bootstrap_api,
	.cuda_device_api_version = nanos6_cuda_device_api,
	.final_api_version = nanos6_final_api,
	.instantiation_api_version = nanos6_instantiation_api,
	.library_mode_api_version = nanos6_library_mode_api,
	.locking_api_version = nanos6_locking_api,
	.polling_api_version = nanos6_polling_api,
	.task_constraints_api_version = nanos6_task_constraints_api,
	.task_execution_api_version = nanos6_task_execution_api,
	.task_info_registration_api_version = nanos6_task_info_registration_api,
	.taskloop_api_version = nanos6_taskloop_api,
	.taskwait_api_version = nanos6_taskwait_api,
	.utils_api_version = nanos6_utils_api
};

void nanos6_preinit(void);

typedef struct {
	pthread_mutex_t _mutex;
	pthread_cond_t _cond;
	int _signaled;
} condition_variable_t;


typedef struct {
	int argc;
	char **argv;
	char **envp;
	int returnCode;
} main_task_args_block_t;

static void main_completion_callback(void *args)
{
	printf("Main finished\n");
	condition_variable_t *condVar = (condition_variable_t *) args;

	pthread_mutex_lock(&condVar->_mutex);
	condVar->_signaled = 1;
	pthread_cond_signal(&condVar->_cond);
	pthread_mutex_unlock(&condVar->_mutex);
}

int start_nanos6(void *start_task, int argc, char **argv)
{
	if (nanos6_check_api_versions(&__user_code_expected_nanos6_api_versions) != 1) {
		return 1;
	}

	// First half of the initialization
	nanos6_preinit();

	condition_variable_t condVar = {
		PTHREAD_MUTEX_INITIALIZER,
		PTHREAD_COND_INITIALIZER,
		0
	};

	// Spawn the main task
	main_task_args_block_t argsBlock = {
		argc,
		argv,
		0
	};

	nanos6_spawn_function(
			start_task,
			&argsBlock,
			main_completion_callback,
			&condVar,
			"start_task");

	// Second half of the initialization
	nanos6_init();

	// Wait for the completion callback
	pthread_mutex_lock(&condVar._mutex);
	while (condVar._signaled == 0)
	{
		printf("Waiting cond wait\n");
		pthread_cond_wait(&condVar._cond, &condVar._mutex);
	}
	pthread_mutex_unlock(&condVar._mutex);

	// Terminate
	nanos6_shutdown();

	return argsBlock.returnCode;
}

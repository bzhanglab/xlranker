<script lang="ts">
  import * as Select from "$lib/components/ui/select/index.js";
  import * as HoverCard from "$lib/components/ui/hover-card/index.js";
  import { Button, buttonVariants } from "$lib/components/ui/button/index.js";
  import FileUpload from "$lib/components/FileUpload.svelte";
  import { AlertCircleIcon, LoaderCircle } from "@lucide/svelte";
  import * as Sheet from "$lib/components/ui/sheet/index.js";
  import * as Alert from "$lib/components/ui/alert/index.js";
  let peptide_file: File[] = $state([]);
  let net_sep = $state("Tab");
  let pep_form: HTMLFormElement;
  let is_loading = $state(false);
  let duplicated = $state(0);
  let err_msg = $state("");

  $effect(() => {
    if (!peptide_file || peptide_file.length === 0) {
      return;
    }
    const file_name = peptide_file ? peptide_file[0].name : "No file selected";
    if (
      file_name.endsWith(".txt") ||
      file_name.endsWith(".tsv") ||
      file_name.endsWith(".csv")
    ) {
      // auto set net_sep based on file extension
      if (file_name.endsWith(".tsv")) {
        net_sep = "Tab";
      } else if (file_name.endsWith(".csv")) {
        net_sep = "Comma";
      }
    } else if (peptide_file && peptide_file.length > 0) {
      // invalid file
      peptide_file = [];
      alert("Please upload a valid .txt, .tsv, or .csv file");
    }
  });

  function upload_network() {
    // post form to /api/upload_network
    is_loading = true;
    duplicated = 0;
    err_msg = "";
    const formData = new FormData(pep_form);
    fetch("/api/upload_network", {
      method: "POST",
      body: formData,
    })
      .then((response) => response.json())
      .then((data) => {
        is_loading = false;
        if (data.status === "ok") {
          console.log("File uploaded successfully:", data);
        } else if (data.status === "err") {
          console.error("Error uploading file:", data.message);
          err_msg = data.message;
        } else if (data.status === "warn") {
          duplicated = data.duplicates;
        }
        // handle success
      })
      .catch((error) => {
        console.error("Error:", error);
        is_loading = false;
        err_msg = error;
        // handle error
      });
  }
</script>

<div class="px-4">
  <h3>Peptide Network</h3>
  Please upload a <Sheet.Root>
    <Sheet.Trigger
      class="decoration-dashed underline hover:decoration-solid font-semibold text-sm"
    >
      Peptide Network
    </Sheet.Trigger>
    <Sheet.Content>
      <Sheet.Header>
        <Sheet.Title>Peptide Network</Sheet.Title>
        <Sheet.Description>
          A text file containing peptide pairs that are observed to be
          cross-linked.
        </Sheet.Description>
      </Sheet.Header>

      <div class="grid flex-1 auto-rows-min gap-6 px-4 text-sm">
        <p>
          The peptide network is the list of peptides that were observed to be
          cross-linked. Each peptide should be the amino acid sequence (e.g. <code
            >KGMFR</code
          >).
        </p>
        <p>
          The file should have two columns, where each column represents a
          peptide. The order of the peptides does not matter and duplicated rows
          are removed automatically.
        </p>
        <p>
          After uploading, XLRanker will verify the format and alert of any
          potential issues.
        </p>

        <h5 class="font-semibold">Example</h5>
        <i>This example is tab-separated</i>
        <pre
          id="example"
          class="whitespace-pre-wrap bg-gray-100 p-4 rounded-md">
EQLAKLMATLR	IIGLDQVAGM
EQLAKLMATLR	KGMFR
EQLAKLMATLR	KGMFR
EQLDNQLDAK	SLNLKHIK
EQLDNQLDAY	SLNLKHIK
</pre>
        <h5 class="font-semibold">Supported Formats</h5>
        <ul class="list-disc list-inside">
          <li>Text files (.txt)</li>
          <li>Comma-separated values (.csv)</li>
          <li>Tab-separated values (.tsv)</li>
        </ul>

        <p>
          Specify the character separating by adjusting the <b
            >Column Separator</b
          > parameter.
        </p>
      </div>
      <Sheet.Footer>
        <Sheet.Close class={buttonVariants({ variant: "outline" })}
          >Close</Sheet.Close
        >
      </Sheet.Footer>
    </Sheet.Content>
  </Sheet.Root> file containing peptide pairs, one pair per line. Each line should
  contain two peptides separated by a tab or comma character.

  <div class="grid gap-2">
    <form
      bind:this={pep_form}
      onsubmit={(e) => {
        e.preventDefault();
        upload_network();
      }}
    >
      <h3>Upload Network</h3>
      <FileUpload
        bind:files={peptide_file}
        name="peptide_network"
        accepted_formats={[".txt", ".csv", ".tsv"]}
      />
      <h4>Column Separator</h4>
      <Select.Root name="col_sep" type="single" bind:value={net_sep}>
        <Select.Trigger class="w-fit min-w-[8rem]">{net_sep}</Select.Trigger>
        <Select.Content>
          <Select.Item value="Tab">Tab</Select.Item>
          <Select.Item value="Comma">Comma</Select.Item>
          <Select.Item value="Space">Space</Select.Item>
        </Select.Content>
      </Select.Root>
      {#if duplicated > 0}
        <Alert.Root class="mt-4 text-yellow-800 bg-yellow-50 border-yellow-200">
          <AlertCircleIcon />
          <Alert.Title>Duplicate Edges Found</Alert.Title>
          <Alert.Description class="mt-2 text-gray-700">
            {`Found and removed ${duplicated} duplicated edge${duplicated > 1 ? "s" : ""} in the network.`}
          </Alert.Description>
        </Alert.Root>
      {/if}
      {#if err_msg !== ""}
        <Alert.Root class="mt-4 bg-red-50 text-red-800 border-red-200">
          <AlertCircleIcon />
          <Alert.Title>Error Uploading File</Alert.Title>
          <Alert.Description class="mt-2 text-gray-700">
            {err_msg}
          </Alert.Description>
        </Alert.Root>
      {/if}
      <div class="flex justify-center">
        <Button type="submit" class="mt-4">
          {#if is_loading}
            <LoaderCircle class="animate-spin" />
            Uploading...
          {:else}
            Submit
          {/if}
        </Button>
      </div>
    </form>
  </div>

  <!-- centered button for submit -->
</div>

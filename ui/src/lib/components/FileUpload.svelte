<script lang="ts">
  import { Label } from "$lib/components/ui/label";
  import { Button } from "$lib/components/ui/button/index.js";
  import * as Card from "$lib/components/ui/card/index.js";

  let { files = $bindable(), multiple_upload = false } = $props();

  let file_text = $state("Click or drag file to this area to upload");

  function handleFileChange(event: Event) {
    const target = event.target as HTMLInputElement;
    if (target.files) {
      files = Array.from(target.files);
    }

    file_text = files.length
      ? files.map((f) => f.name).join(", ")
      : "Click or drag file to this area to upload";
  }

  function handleDrop(event: DragEvent) {
    event.preventDefault();
    if (event.dataTransfer?.files) {
      files = Array.from(event.dataTransfer.files);
    }
  }

  function handleDragOver(event: DragEvent) {
    event.preventDefault();
  }
</script>

<!-- Card container for visual structure -->
<div class="flex justify-center">
  <div class="w-full max-w-md">
    <!-- Drag-and-drop upload area -->
    <div
      class="content-center mx-auto border-1 text-gray-600 border-dashed border-gray-300 rounded-md p-6 text-center cursor-pointer hover:border-gray-500 hover:text-black transition"
      ondrop={handleDrop}
      ondragover={handleDragOver}
      role="button"
      tabindex="0"
    >
      <Label for="file-upload" class="cursor-pointer flex flex-col">
        <p class="text-md text-center">{file_text}</p>
      </Label>
      <input
        id="file-upload"
        type="file"
        class="hidden"
        onchange={handleFileChange}
        multiple={multiple_upload}
        bind:value={files}
      />
    </div>
  </div>
</div>
